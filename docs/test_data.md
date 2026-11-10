# Test run on a small short-read sample

This test runs the short-read part of HitFlow on public data: QC, alignment,
BQSR, GVCF calling, single-sample filtering, Manta and smoove SV calling with
SURVIVOR merging, and MultiQC. It checks that the software environments and
cluster setup work before you download the full reference data for a real run.

It needs about 100 MB of downloads and no annotation databases.

## Test data

`scripts/pipeline/prepare_test_data.sh` reads only the needed regions from
public indexed files:

- `data/HG002_chr20_R{1,2}.fq.gz`: about 70,000 read pairs (2 × 250 bp) of the
  Genome in a Bottle sample HG002 (NA24385), from the NIST Illumina 2 × 250 bp
  data set aligned to GRCh38 by GIAB (`HG002.GRCh38.2x250.bam` in
  <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/>).
  These are the pairs whose reads both have their primary alignment in
  chr20:63,000,000–63,500,000 (about 30× coverage). The region contains 68 SVs
  of at least 50 bp in the HG002 Q100 v1.1 benchmark, so the SV callers have
  real events to find. HitFlow aligns the reads again with BWA-MEM.
- `ref/chr20.fa`: chr20 of `GRCh38_full_analysis_set_plus_decoy_hla.fa` (1000
  Genomes FTP), the reference used by the full configuration. The workflow
  builds the FASTA index, sequence dictionary and BWA index.
- `ref/known_sites.chr20.vcf.gz`: known variant sites for BQSR in
  chr20:62.9–63.6 Mb (about 24,000 sites) from the 1000 Genomes high-coverage
  GRCh38 call set.

Because the reads cover only half a megabase and no databases are installed:

- Per-chromosome steps run on chr20 only.
- Small-variant annotation (`ngs_single_sample_variants.annotation.tool: "none"`),
  AnnotSV and VerifyBamID2 are switched off, and smoove runs without an exclude
  BED.
- The genome-wide coverage check (Qualimap) fails, and the sample is listed in
  `results/ngs/qc_reports/failed_samples.txt`. This is expected. A Qualimap
  failure is only recorded, so the analysis continues.

The test checks that every step runs; it does not measure accuracy. For
accuracy, see the HG002 configurations in [config/examples/](../config/examples/).

## Requirements

- Everything in the README's [Requirements](../README.md#requirements).
- For `prepare_test_data.sh`: internet access, and `samtools`, `tabix` and
  `bgzip` from htslib 1.10 or later with network support. The samtools
  container that HitFlow uses has all three (see below).

## Run

Use a separate run directory:

```bash
HITFLOW=/path/to/HitFlow
mkdir hitflow_test && cd hitflow_test
cp $HITFLOW/docs/test_data/config.yaml $HITFLOW/docs/test_data/samples_ngs.csv .

# Writes data/ and ref/ in the current directory (about 5 minutes).
# With SingularityCE, replace "apptainer exec" with "singularity exec".
apptainer exec -B "$HITFLOW" docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0 \
    bash $HITFLOW/scripts/pipeline/prepare_test_data.sh
# or, with samtools, tabix and bgzip installed:
# bash $HITFLOW/scripts/pipeline/prepare_test_data.sh

snakemake -s $HITFLOW/Snakefile -n    # dry run: 45 jobs
snakemake -s $HITFLOW/Snakefile --profile $HITFLOW/profiles/local-apptainer --cores 8 \
    --apptainer-args "-B $HITFLOW --cleanenv"
# or on Slurm:
snakemake -s $HITFLOW/Snakefile --profile $HITFLOW/profiles/slurm-apptainer \
    --apptainer-args "-B $HITFLOW --cleanenv"
```

If the internet is reached through a proxy, set `https_proxy` and `http_proxy`
before running the script.

## Expected output

- `results/multiqc/multiqc_report.html`
- `results/ngs/alignment/HG002/bam/HG002.dup.sorted.recal.bam`
- `results/ngs/variants/single_sample/HG002/variant_calling/normalized/HG002.filtered.normalized.vcf.gz`
- `results/ngs/structural_variants/single_sample/HG002/merged/HG002.ngs_sv.merged.vcf.gz`

Our run on a Slurm cluster (45 jobs, about 9 minutes including queue time)
gave 1,340 filtered small variants, 9 Manta SVs, 3 smoove SVs and 3 merged
SVs, with HG002 listed in `failed_samples.txt` by Qualimap as described above.
Your counts may differ slightly with other software versions.
