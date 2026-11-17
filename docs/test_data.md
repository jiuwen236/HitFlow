# Test run on a small short-read and HiFi sample

This test runs HitFlow on public data: short-read QC, alignment, BQSR, GVCF
calling and single-sample filtering; Manta and smoove SV calling with SURVIVOR
merging; HiFi alignment and SV calling with Sniffles2, cuteSV, SVIM and pbsv;
HitSV on the short reads, on the HiFi reads and on both together (Hybrid);
callset-level integration of the short-read and HiFi call sets; and MultiQC. It checks that the software environments and cluster setup work before
you download the full reference data for a real run.

It needs about 200 MB of downloads and no annotation databases.

## Test data

`scripts/pipeline/prepare_test_data.sh` reads only the needed regions from
public indexed files:

All reads come from chr20:61,000,000 to the end of the chromosome
(64,444,167), the last 3.44 Mb. This region holds 326 SVs of at least 50 bp
(123 deletions, 203 insertions) in the HG002 Q100 v1.1 benchmark, so every SV
caller has real events to find.

- `data/HG002_chr20_R{1,2}.fq.gz`: about 465,000 read pairs (2 × 250 bp) of the
  Genome in a Bottle sample HG002 (NA24385), from the NIST Illumina 2 × 250 bp
  data set aligned to GRCh38 by GIAB (`HG002.GRCh38.2x250.bam` in
  <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/>).
  These are the pairs whose reads both have their primary alignment in the
  region. HitFlow aligns the reads again with BWA-MEM.
- `data/HG002_chr20_hifi.fq.gz`: about 780 PacBio HiFi (Revio) reads of the same
  sample and region, from `HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam`
  in <https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/>,
  subsampled to about 4×. HitSV Hybrid combines deep short-read data with
  shallow long-read data, which is what this pair of inputs is for. HitFlow
  aligns the reads again with minimap2. `pairs.csv` links this sample to the
  short-read one, which both HitSV Hybrid and the callset-level integration
  (`enable_cross_platform_sv`) need.
- `ref/chr20.fa`: chr20 of `GRCh38_full_analysis_set_plus_decoy_hla.fa` (1000
  Genomes FTP), the reference used by the full configuration. The workflow
  builds the FASTA index, sequence dictionary and BWA index.
- `ref/known_sites.chr20.vcf.gz`: known variant sites for BQSR in
  chr20:60.9 Mb–end (about 112,000 sites) from the 1000 Genomes high-coverage
  GRCh38 call set.

Because the reads cover only the last 3.44 Mb of chr20 and no databases are installed:

- Per-chromosome steps run on chr20 only.
- Small-variant annotation (`ngs_single_sample_variants.annotation.tool: "none"`),
  AnnotSV, VerifyBamID2, CNVpytor and NanoVar (ONT only) are switched off, and
  smoove runs without an exclude BED.
- The genome-wide coverage check (Qualimap) fails, and the sample is listed in
  `results/ngs/qc_reports/failed_samples.txt`. This is expected. A Qualimap
  failure is only recorded, so the analysis continues.
- HitSV runs with `repair_padding: true` and
  `invalid_allele_policy: "quarantine"`, as in
  [config/examples/hg002_hybrid_sv/](../config/examples/hg002_hybrid_sv/). The
  repairs and any set-aside alleles are counted in
  `results/hitsv_work/*/*/vcf_audit.json`.
- MultiQC plot export (`multiqc.export_plots`) is off because it can hang on a
  compute node; the HTML report still has the interactive plots.

## Requirements

- Everything in the README's [Requirements](../README.md#requirements).
- For `prepare_test_data.sh`: internet access, and `samtools`, `tabix` and
  `bgzip` from htslib 1.10 or later with network support. The samtools
  container that HitFlow uses has all three (see below).
- The HitSV binary in `bin/`, installed by
  `scripts/pipeline/install_hitsv.sh`. To run the test without it, set the
  three `hitsv` switches in `config.yaml` (`ngs_sv.hitsv.enabled`,
  `tgs_config.sv_detection.hitsv.enabled` and `hitsv.hybrid.enabled`) to
  `false`.
- Manta is downloaded from GitHub into `bin/manta/` by the first Manta job if it
  is not there yet. Without direct GitHub access from that job, install it
  beforehand with `bash -c "source scripts/pipeline/import_manta.sh"` from the
  HitFlow directory.

## Run

Use a separate run directory:

```bash
HITFLOW=/path/to/HitFlow
mkdir hitflow_test && cd hitflow_test
cp $HITFLOW/docs/test_data/config.yaml $HITFLOW/docs/test_data/samples_ngs.csv \
   $HITFLOW/docs/test_data/samples_tgs.csv $HITFLOW/docs/test_data/pairs.csv .

# Writes data/ and ref/ in the current directory (about 15 minutes).
# With SingularityCE, replace "apptainer exec" with "singularity exec".
apptainer exec -B "$HITFLOW" docker://quay.io/biocontainers/samtools:1.20--h50ea8bc_0 \
    bash $HITFLOW/scripts/pipeline/prepare_test_data.sh
# or, with samtools, tabix and bgzip installed:
# bash $HITFLOW/scripts/pipeline/prepare_test_data.sh

snakemake -s $HITFLOW/Snakefile -n    # dry run: 92 jobs
snakemake -s $HITFLOW/Snakefile --profile $HITFLOW/profiles/local-apptainer --cores 8 \
    --apptainer-args "-B $HITFLOW --cleanenv"
# or on Slurm:
snakemake -s $HITFLOW/Snakefile --profile $HITFLOW/profiles/slurm-apptainer \
    --apptainer-args "-B $HITFLOW --cleanenv"
```

If the internet is reached through a proxy, set `https_proxy` and `http_proxy`
before running the script.

With Apptainer, run Snakemake from the real path of the run directory
(`cd "$(pwd -P)"`) if the path contains a symbolic link, for example a home
directory that links to another disk. Otherwise the jobs cannot find their
relative paths inside the container.

## Expected output

- `results/multiqc/multiqc_report.html`
- `results/ngs/alignment/HG002/bam/HG002.dup.sorted.recal.bam`
- `results/ngs/variants/single_sample/HG002/variant_calling/normalized/HG002.filtered.normalized.vcf.gz`
- `results/ngs/structural_variants/single_sample/HG002/merged/HG002.ngs_sv.merged.vcf.gz`
- `results/tgs/structural_variants/single_sample/HG002_HiFi/final/HG002_HiFi.sv.final.vcf.gz`
- `results/cross_platform/hitsv/HG002/HG002.hitsv.vcf.gz`
- `results/cross_platform/structural_variants/HG002/HG002.cross_sv.merged.vcf.gz`
  and `HG002.concordance.tsv`

Our runs on a Slurm cluster (CentOS 7) and on one Debian 12 machine gave the
same 7,562 filtered small variants and these SV counts:

| Callset | Input | SVs |
| --- | --- | --- |
| Manta | short reads | 29 |
| smoove | short reads | 11 |
| SURVIVOR merge of both | short reads | 25 |
| HitSV | short reads | 49 |
| Sniffles2, cuteSV, SVIM and pbsv merged | HiFi | 185 |
| HitSV | HiFi | 245 |
| HitSV Hybrid | both | 218 |
| SURVIVOR merge of both platforms (callset level) | both | 149 |

The callset-level merge keeps every site, so its concordance table
(`HG002.concordance.tsv`) splits those 149 SVs into 20 found on both platforms,
5 short-read only and 124 HiFi only.

HG002 is listed in `failed_samples.txt` by Qualimap as described above. Your
counts may differ slightly with other software versions. The test checks that
every step runs; for accuracy numbers, see the HG002 configurations in
[config/examples/](../config/examples/).

## Resources

Measured on a Debian 12 machine with 8 cores and 16 GB of memory, with
`profiles/local-apptainer` and `--cores 8`:

| Step | Time | Download | Disk |
| --- | --- | --- | --- |
| Miniforge and the `snakemake` environment | 5 minutes | – | 2.2 GB |
| `prepare_test_data.sh` | 10 minutes | 0.2 GB | 0.2 GB |
| Pulling 24 container images (first run only) | 33 minutes | 5.2 GB | 4.6 GB, plus 9.5 GB in `~/.apptainer/cache` |
| The 92 jobs | 10 minutes | – | 0.5 GB in `results/`, 1.1 GB in `tmp/` |

The jobs used about half a CPU-hour. The largest job was GATK HaplotypeCaller
(2.8 GB, 3 minutes), and the whole run never used more than 6 GB of memory.
`apptainer cache clean` frees the image cache once the images are in
`.snakemake/singularity/`.
