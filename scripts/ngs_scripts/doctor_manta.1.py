#!/usr/bin/env python3
# Source: adapted from the Hall lab sv-pipeline (https://github.com/hall-lab/sv-pipeline).
"""
This script normalizes Manta VCF records into a form compatible with downstream
SV merging/annotation.

Inputs
  - stdin or -i/--vcf <path>
  - -m/--max_ins <int>: maximum insertion length to cap SVLEN for INS (default 1000)

Outputs
  - VCF to stdout

Notes
  - Only single-sample VCFs are expected (as produced by Manta per-sample).
  - Header lines for PRPOS, PREND, STRANDS, SU, PE, SR, INSLEN_ORIG, CIPOS95,
    CIEND95, SECONDARY are added if they do not already exist.
"""

import sys
import argparse
import re


def add_arguments_to_parser(parser):
    parser.add_argument('-i', '--vcf', metavar='<VCF>', dest='manta_vcf', help='manta input vcf')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true')
    parser.add_argument('-s', '--slop', dest='slop', default=0, required=False, help='padding to either side')
    parser.add_argument('-m', '--max_ins', dest='max_ins', default=1000, type=int, required=False, help='maximum insert size')


def command_parser():
    parser = argparse.ArgumentParser(description='Normalize Manta VCF to downstream-friendly representation')
    add_arguments_to_parser(parser)
    return parser


def parse_info(info_str):
    if info_str == '.' or not info_str:
        return {}
    d = {}
    for ent in info_str.split(';'):
        if not ent:
            continue
        if '=' in ent:
            k, v = ent.split('=', 1)
            d[k] = v
        else:
            d[ent] = True
    return d


def info_to_str(info):
    if not info:
        return '.'
    parts = []
    for k, v in info.items():
        if v is True:
            parts.append(k)
        else:
            parts.append(f'{k}={v}')
    return ';'.join(parts)


def parse_format(fmt, sample):
    keys = fmt.split(':') if fmt else []
    vals = sample.split(':') if sample else []
    return dict(zip(keys, vals))


def get_pr_sr(fmt_map):
    pe = 0
    sr = 0
    if 'PR' in fmt_map:
        try:
            pe = int(fmt_map['PR'].split(',')[1])
        except Exception:
            pass
    if 'SR' in fmt_map:
        try:
            sr = int(fmt_map['SR'].split(',')[1])
        except Exception:
            pass
    return pe, sr


def parse_ci(ci_str):
    try:
        a, b = ci_str.split(',')
        return int(a), int(b)
    except Exception:
        return 0, 0


def uniform_pr_str(start, stop):
    length = stop - start + 1
    if length <= 0:
        return '1.0'
    val = 1.0 / float(length)
    return ','.join(str(val) for _ in range(length))


def parse_bnd_partner(alt):
    m = re.search(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]', alt)
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def convert_variant(cols, max_ins):
    chrom = cols[0]
    pos = int(cols[1])
    ref = cols[3]
    alt = cols[4]
    info = parse_info(cols[7])
    fmt = cols[8] if len(cols) > 8 else None
    sample = cols[9] if len(cols) > 9 else None

    # Read counts from FORMAT -> PE, SR, SU
    pe = sr = 0
    if fmt and sample:
        fmap = parse_format(fmt, sample)
        pe, sr = get_pr_sr(fmap)
        info['PE'] = str(pe)
        info['SR'] = str(sr)
        info['SU'] = str(pe + sr)

    # CIPOS/CIEND -> PRPOS/PREND; duplicate into CIPOS95/CIEND95
    imprec = False
    if 'CIPOS' in info:
        start, stop = parse_ci(info['CIPOS'])
        info['PRPOS'] = uniform_pr_str(start, stop)
        imprec = True
    else:
        info['CIPOS'] = '0,0'
        info['PRPOS'] = '1.0'
    if 'CIEND' in info:
        start, stop = parse_ci(info['CIEND'])
        info['PREND'] = uniform_pr_str(start, stop)
        imprec = True
    else:
        info['CIEND'] = '0,0'
        info['PREND'] = '1.0'
    info['CIPOS95'] = info.get('CIPOS', '0,0')
    info['CIEND95'] = info.get('CIEND', '0,0')
    if imprec:
        info['IMPRECISE'] = True

    su = info.get('SU')  # may be string or missing
    svtype = info.get('SVTYPE')

    if svtype == 'DEL':
        ref = 'N'
        alt = '<DEL>'
        if su is not None:
            info['STRANDS'] = f'+-:{su}'
    elif svtype == 'DUP':
        ref = 'N'
        alt = '<DUP>'
        if su is not None:
            info['STRANDS'] = f'-+:{su}'
    elif svtype == 'INV':
        ref = 'N'
        end = info.get('END')
        if end is not None:
            if 'INV3' in info:
                alt = f'N]{chrom}:{end}]'
                if su is not None:
                    info['STRANDS'] = f'++:{su}'
            else:
                # Match original behavior exactly
                alt = f'[{chrom}:{end}['
                if su is not None:
                    info['STRANDS'] = f'--:{su}'
            info['SVTYPE'] = 'BND'
    elif svtype == 'INS':
        ref = 'N'
        alt = '<INS>'
        if su is not None:
            info['STRANDS'] = f'+.:{su}'
        orig_len = '.'
        new_len = max_ins
        if 'SVLEN' in info:
            try:
                svlen = int(info['SVLEN'])
                orig_len = str(svlen)
                if svlen < max_ins:
                    new_len = svlen
            except Exception:
                pass
        info['SVLEN'] = str(new_len)
        info['INSLEN_ORIG'] = str(orig_len)
    elif svtype == 'BND':
        ref = 'N'
        ff = alt.find('[')
        if ff == 0:
            strands = '--:'
            ff1 = alt.find('[', 1)
            newalt = alt[0:(ff1 + 1)] + 'N' if ff1 != -1 else alt + 'N'
        elif ff > 0:
            strands = '+-:'
            newalt = 'N' + alt[ff:]
        else:
            ff = alt.find(']')
            if ff == 0:
                strands = '-+:'
                ff1 = alt.find(']', 1)
                newalt = alt[0:(ff1 + 1)] + 'N' if ff1 != -1 else alt + 'N'
            else:
                strands = '++:'
                newalt = 'N' + alt[ff:]
        alt = newalt
        if su is not None:
            info['STRANDS'] = strands + str(su)
        chrom2, pos2 = parse_bnd_partner(alt)
        if chrom2 is not None and pos2 is not None:
            is_secondary = (chrom2 < chrom) or (chrom2 == chrom and int(pos2) < pos)
            if is_secondary:
                info['SECONDARY'] = True

    cols[3] = ref
    cols[4] = alt
    cols[7] = info_to_str(info)
    return cols


def header_has_info(header_lines, info_id):
    pat = re.compile(r'^##INFO=<ID=' + re.escape(info_id) + r'[>,]')
    return any(pat.search(h) for h in header_lines)


def emit_augmented_header(header_lines, out):
    additions = [
        '##INFO=<ID=PRPOS,Number=1,Type=String,Description="Breakpoint probability dist">',
        '##INFO=<ID=PREND,Number=1,Type=String,Description="Breakpoint probability dist">',
        '##INFO=<ID=STRANDS,Number=.,Type=String,Description="Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--">',
        '##INFO=<ID=SU,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">',
        '##INFO=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">',
        '##INFO=<ID=SR,Number=.,Type=Integer,Description="Number of split reads supporting the variant across all samples">',
        '##INFO=<ID=INSLEN_ORIG,Number=.,Type=Integer,Description="Original insertion length">',
        '##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description="Confidence interval (95%) around POS for imprecise variants">',
        '##INFO=<ID=CIEND95,Number=2,Type=Integer,Description="Confidence interval (95%) around END for imprecise variants">',
        '##INFO=<ID=SECONDARY,Number=0,Type=Flag,Description="Secondary breakend in a multi-line variant">',
    ]
    for add in additions:
        info_id = add.split('ID=')[1].split(',')[0].split('>')[0]
        if not header_has_info(header_lines, info_id):
            header_lines.append(add + '\n')
    for h in header_lines:
        out.write(h)


def run_from_args(args):
    infile = open(args.manta_vcf, 'rt') if args.manta_vcf else sys.stdin
    out = sys.stdout

    header_lines = []
    for line in infile:
        if line.startswith('##'):
            header_lines.append(line)
            continue
        if line.startswith('#CHROM'):
            # Inject additional INFO definitions before writing header
            emit_augmented_header(header_lines, out)
            out.write(line)
            break
    else:
        # No #CHROM line found; just emit collected header and return
        emit_augmented_header(header_lines, out)
        return

    # Process body
    for line in infile:
        if not line or line[0] == '#':
            out.write(line)
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 8:
            # Not a valid VCF row; emit unchanged
            out.write(line)
            continue
        cols = convert_variant(cols, args.max_ins)
        out.write('\t'.join(cols) + '\n')


if __name__ == '__main__':
    parser = command_parser()
    args = parser.parse_args()
    run_from_args(args)
