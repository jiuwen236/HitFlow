#!/usr/bin/env python3
import sys
import os
import argparse
import math
import gzip
from datetime import datetime


def open_text(path):
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def write_header(out, reference, sample):
    out.write('##fileformat=VCFv4.2\n')
    out.write('##fileDate=%s\n' % datetime.now().strftime('%Y%m%d'))
    if reference:
        out.write('##reference=%s\n' % reference)
    out.write('##source=CNVpytor\n')
    out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
    out.write('##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n')
    out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
    out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    out.write('##INFO=<ID=natorRD,Number=1,Type=Float,Description="Normalized RD">\n')
    out.write('##INFO=<ID=natorP1,Number=1,Type=Float,Description="p-val by t-test">\n')
    out.write('##INFO=<ID=natorP2,Number=1,Type=Float,Description="p-val by Gaussian tail">\n')
    out.write('##INFO=<ID=natorP3,Number=1,Type=Float,Description="p-val by t-test (middle)">\n')
    out.write('##INFO=<ID=natorP4,Number=1,Type=Float,Description="p-val by Gaussian tail (middle)">\n')
    out.write('##INFO=<ID=natorQ0,Number=1,Type=Float,Description="Fraction of reads with 0 mapping quality">\n')
    out.write('##ALT=<ID=DEL,Description="Deletion">\n')
    out.write('##ALT=<ID=DUP,Description="Duplication">\n')
    out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')


def parse_summary_line(line):
    parts = line.strip().split()
    if len(parts) < 9:
        return None
    svtype = parts[0].lower()
    coord = parts[1]
    try:
        chrom, rest = coord.split(':', 1)
        start_str, end_str = rest.split('-', 1)
        start = int(start_str)
        end = int(end_str)
    except Exception:
        return None

    # Expected CNVpytor summary columns: type, chr:start-end, length, RD, p1, p2, p3, p4, q0 [, pe]
    length_val = parts[2]
    rd = parts[3]
    p1 = parts[4] if len(parts) > 4 else ''
    p2 = parts[5] if len(parts) > 5 else ''
    p3 = parts[6] if len(parts) > 6 else ''
    p4 = parts[7] if len(parts) > 7 else ''
    q0 = parts[8] if len(parts) > 8 else ''

    is_del = (svtype == 'deletion')
    is_dup = (svtype == 'duplication')
    if not (is_del or is_dup):
        return None

    try:
        svlen = int(length_val)
    except Exception:
        svlen = abs(end - start)
    if is_del:
        svlen = -abs(svlen)
    else:
        svlen = abs(svlen)

    info_fields = [f'END={end}']
    info_fields.append('SVTYPE=%s' % ('DEL' if is_del else 'DUP'))
    info_fields.append('SVLEN=%d' % svlen)
    info_fields.append('IMPRECISE')
    def _append_float(tag, val):
        try:
            if val != '' and val.lower() != 'nan':
                float(val)
                info_fields.append('%s=%s' % (tag, val))
        except Exception:
            pass
    _append_float('natorRD', rd)
    _append_float('natorP1', p1)
    _append_float('natorP2', p2)
    _append_float('natorP3', p3)
    _append_float('natorP4', p4)
    _append_float('natorQ0', q0)

    alt = '<DEL>' if is_del else '<DUP>'

    return chrom, start, alt, ';'.join(info_fields)


def convert(summary_path, sample, reference, ref_base):
    out = sys.stdout
    write_header(out, reference, sample)
    count = 0
    with open_text(summary_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parsed = parse_summary_line(line)
            if parsed is None:
                continue
            chrom, pos, alt, info = parsed
            count += 1
            vid = f'CNV_{count}'
            out.write(f'{chrom}\t{pos}\t{vid}\t{ref_base}\t{alt}\t.\tPASS\t{info}\n')


def main():
    p = argparse.ArgumentParser(description='Convert CNVpytor summary to VCF (stdout).')
    p.add_argument('--sample', required=True, help='Sample ID (for downstream use; not embedded here)')
    p.add_argument('--input', required=True, help='CNVpytor summary file (text, possibly gz)')
    p.add_argument('--reference', default='', help='Reference name/path to record in VCF header')
    p.add_argument('--ref-base', default='.', help='REF base to emit in VCF (default ".")')
    args = p.parse_args()

    if not os.path.exists(args.input):
        sys.stderr.write(f'Input not found: {args.input}\n')
        return 2

    ref_base = args.ref_base
    if not ref_base:
        ref_base = '.'
    if len(ref_base) > 1 and ref_base != '.':
        # Keep simple and robust: enforce single base or '.'
        ref_base = ref_base[0]
    convert(args.input, args.sample, args.reference, ref_base)


if __name__ == '__main__':
    sys.exit(main() or 0)


