#!/usr/bin/env bash
# Download the small public data set for the test run (docs/test_data.md):
#   data/HG002_chr20_R{1,2}.fq.gz   GIAB HG002 Illumina 2x250 bp read pairs aligned
#                                   to GRCh38 chr20:61.0 Mb-end
#   data/HG002_chr20_hifi.fq.gz     GIAB HG002 PacBio HiFi (Revio) reads from the
#                                   same region, subsampled to about 4x
#   ref/chr20.fa                    chr20 of the GRCh38 analysis set used by the
#                                   full configuration
#   ref/known_sites.chr20.vcf.gz    known variant sites for BQSR around the test
#                                   region (1000 Genomes high-coverage call set)
# Only the needed regions are read from the remote indexed files (about 400 MB in
# total). The FASTA index, sequence dictionary and BWA index are built by the
# workflow.
#
# Usage: prepare_test_data.sh [DEST]     (default DEST: current directory)
# Requires samtools, tabix and bgzip (htslib >= 1.10, with network support), e.g.
# the samtools container used by HitFlow:
#   apptainer exec docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0 \
#       bash prepare_test_data.sh
# Set https_proxy/http_proxy if the internet is reached through a proxy.
set -euo pipefail

dest=$(cd "${1:-.}" && pwd)
mkdir -p "$dest/ref" "$dest/data"

BAM_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam
FASTA_URL=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
SITES_URL=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr20.filtered.shapeit2-duohmm-phased.vcf.gz
HIFI_URL=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
CHR20_LENGTH=64444167
READS_REGION=chr20:61000000-$CHR20_LENGTH
SITES_REGION=chr20:60900000-$CHR20_LENGTH
# samtools view -s SEED.FRACTION: keep about 4x of the HiFi data (about 34x
# in this region).
HIFI_SUBSAMPLE=11.12

for tool in samtools tabix bgzip; do
    command -v "$tool" >/dev/null || { echo "Error: $tool is required" >&2; exit 1; }
done

# Remote index files may be cached in the working directory, so work in a temp dir.
tmp=$(mktemp -d "$dest/.prepare_test_data.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

r1="$dest/data/HG002_chr20_R1.fq.gz"
r2="$dest/data/HG002_chr20_R2.fq.gz"
if [[ -s "$r1" && -s "$r2" ]]; then
    echo "Found $r1 and $r2"
else
    echo "Fetching HG002 read pairs in $READS_REGION (a few minutes) ..."
    # Primary alignments only, so each read appears once; pairs with a mate
    # outside the region are dropped by samtools fastq.
    (cd "$tmp" && samtools view -b -F 0x900 -o region.bam "$BAM_URL" "$READS_REGION")
    samtools sort -n -o "$tmp/region.qsort.bam" "$tmp/region.bam"
    samtools fastq -n -F 0x900 -1 "$tmp/R1.fq.gz" -2 "$tmp/R2.fq.gz" \
        -0 /dev/null -s /dev/null "$tmp/region.qsort.bam"
    mv "$tmp/R1.fq.gz" "$r1"
    mv "$tmp/R2.fq.gz" "$r2"
fi

hifi="$dest/data/HG002_chr20_hifi.fq.gz"
if [[ -s "$hifi" ]]; then
    echo "Found $hifi"
else
    echo "Fetching HG002 HiFi reads in $READS_REGION (a few minutes) ..."
    (cd "$tmp" && samtools view -b -F 0x900 -s "$HIFI_SUBSAMPLE" -o hifi.bam "$HIFI_URL" "$READS_REGION")
    samtools fastq -F 0x900 "$tmp/hifi.bam" | bgzip -c > "$tmp/hifi.fq.gz"
    mv "$tmp/hifi.fq.gz" "$hifi"
fi

if [[ -s "$dest/ref/chr20.fa" ]]; then
    echo "Found $dest/ref/chr20.fa"
else
    echo "Fetching GRCh38 chr20 ..."
    (cd "$tmp" && samtools faidx "$FASTA_URL" chr20 > chr20.fa)
    mv "$tmp/chr20.fa" "$dest/ref/chr20.fa"
fi

sites="$dest/ref/known_sites.chr20.vcf.gz"
if [[ -s "$sites" && -s "$sites.tbi" ]]; then
    echo "Found $sites"
else
    echo "Fetching known sites in $SITES_REGION ..."
    # Sites only (first 8 columns); INFO is dropped because BQSR uses positions only.
    (cd "$tmp" && tabix "$SITES_URL" "$SITES_REGION") \
        | awk 'BEGIN { OFS = "\t" } { print $1, $2, $3, $4, $5, ".", $7, "." }' > "$tmp/body.txt"
    [[ -s "$tmp/body.txt" ]] || { echo "Error: no sites returned from $SITES_URL" >&2; exit 1; }
    {
        printf '##fileformat=VCFv4.2\n'
        printf '##contig=<ID=chr20,length=%s>\n' "$CHR20_LENGTH"
        printf '##source=1000 Genomes high-coverage GRCh38 (%s)\n' "$SITES_REGION"
        printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        cat "$tmp/body.txt"
    } | bgzip -c > "$tmp/sites.vcf.gz"
    mv "$tmp/sites.vcf.gz" "$sites"
    tabix -f -p vcf "$sites"
    echo "Wrote $(wc -l < "$tmp/body.txt") sites to $sites"
fi
echo "Test data ready in $dest/data and $dest/ref"
