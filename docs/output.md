# Output files

Paths below use the default output settings (`output.root: results`, logs in
`logs/`, temporary files in `tmp/`). `{sample}` is the sample name and
`{chrom}` a chromosome. Which files are produced depends on the `workflows.*`
switches, the tool switches and the input type. All output directories can be
changed under `output` in the configuration.

## Reports

| File | Content |
| --- | --- |
| `results/multiqc/multiqc_report.html` | MultiQC report (fastp, FastQC, Qualimap, VerifyBamID2, bcftools stats, snpEff) |
| `results/multiqc/multiqc_report_data/`, `multiqc_report_plots/` | Parsed data and exported plots |
| `results/multiqc/multiqc_inputs.txt` | Files that were passed to MultiQC |
| `results/env_check/report.tsv` | Environment check summary (`config/env_check.yaml`); per tool in `results/env_check/tools/{tool}.tsv` |

The MultiQC version used (1.33) has no NanoPlot module, so long-read QC is in
the per-sample NanoPlot reports.

## Short reads: preprocessing and QC

| File | Content |
| --- | --- |
| `results/ngs/preprocessing/{sample}/fastp/{sample}.clean.R{1,2}.fq.gz` | Trimmed reads |
| `results/ngs/preprocessing/{sample}/fastp/{sample}.fastp.{html,json}` | fastp report |
| `results/ngs/preprocessing/{sample}/fastqc/{sample}_{1,2}_fastqc.{html,zip}` | FastQC reports |
| `results/ngs/preprocessing/{sample}/fastqc/{sample}.fastqc.report` | FastQC pass/fail summary used by the QC filter |
| `results/ngs/alignment/{sample}/qc/qualimap/qualimapReport.html` | Qualimap report; summary in `genome_results.txt` |
| `results/ngs/alignment/{sample}/qc/verifybamid/{sample}.selfSM` | VerifyBamID2 contamination estimate |
| `results/ngs/qc_reports/` | Copies of the QC reports of samples that passed |
| `results/ngs/qc_reports/failed_samples.txt` | Samples that failed a QC threshold; their reports are in `failed_samples/` |

QC thresholds are set in `scripts/ngs_scripts/configure.v1` (WGS) or
`configure.wes` (WES), selected by `configure_path`. A sample that fails the
FastQC or VerifyBamID2 check stops at that step (the job fails; other samples
continue with `--keep-going`). A sample that fails the Qualimap check is only
listed, and its analysis continues.

## Short reads: alignment

| File | Content |
| --- | --- |
| `results/ngs/alignment/{sample}/bam/{sample}.dup.sorted.bam` | Sorted BAM with duplicates marked |
| `results/ngs/alignment/{sample}/bam/{sample}.dup.sorted.recal.bam` | BAM after BQSR (input of variant calling) |
| `results/ngs/alignment/{sample}/bam/{sample}.{split,discordant}.bam` | Split and discordant reads (samblaster) |
| `results/ngs/alignment/{sample}/bqsr/{sample}.dup.sorted.recalibration_report.grp` | BQSR table |

## Short reads: SNVs and indels

| File | Content |
| --- | --- |
| `results/ngs/variants/single_sample/{sample}/gvcf/{sample}.{chrom}.hc.g.vcf.bgz` | Per-chromosome GVCFs |
| `results/ngs/variants/single_sample/{sample}/variant_calling/normalized/{sample}.filtered.normalized.vcf.gz` | Single-sample hard-filtered, normalised VCF |
| `results/ngs/variants/single_sample/{sample}/variant_calling/annotated/{sample}.funcotator.{vcf.gz,tsv}` | Funcotator annotation |
| `results/ngs/variants/single_sample/{sample}/variant_calling/annotated/{sample}.snpeff.{vcf.gz,tsv}` | snpEff annotation; HTML summary in `{sample}.snpeff_summary.html` |
| `results/ngs/variants/cohort/genotyped/{chrom}.genotyped.vcf.gz` | Joint-genotyped cohort VCFs |
| `results/ngs/variants/cohort/vqsr/cohort.vqsr.vcf.gz` | Cohort VCF after VQSR |
| `results/ngs/variants/cohort/filtered/cohort.filtered.recode.vcf.gz` | Filtered cohort VCF |
| `results/ngs/variants/cohort/phased/{chrom}.phased.vcf.gz` | SHAPEIT5-phased cohort VCFs |

## Short reads: structural variants

| File | Content |
| --- | --- |
| `results/ngs/structural_variants/single_sample/{sample}/manta/{sample}.manta.vcf.gz` | Manta calls; `{sample}.manta.doctored.vcf.gz` is the normalised version |
| `results/ngs/structural_variants/single_sample/{sample}/smoove/{sample}-smoove.genotyped.vcf.gz` | smoove/LUMPY calls |
| `results/ngs/structural_variants/single_sample/{sample}/cnvpytor/{sample}.cnvpytor.vcf.gz` | CNVpytor calls (optional) |
| `results/ngs/structural_variants/single_sample/{sample}/hitsv/{sample}.hitsv.vcf.gz` | HitSV short-read calls (optional) |
| `results/ngs/structural_variants/single_sample/{sample}/merged/{sample}.ngs_sv.merged.vcf.gz` | Per-sample merged SVs |
| `results/ngs/structural_variants/single_sample/{sample}/annotation/{sample}.ngs.annotsv.tsv` | AnnotSV annotation of the merged SVs |
| `results/ngs/structural_variants/cohort/merged/cohort.ngs_sv.merged.vcf.gz` | Cohort SV merge |

## Long reads: alignment and QC

| File | Content |
| --- | --- |
| `results/tgs/preprocessing/{sample}/alignment/{sample}.{platform}.sorted.bam` | Sorted BAM; `{platform}` is `ont` or `pb` |
| `results/tgs/preprocessing/{sample}/qc/nanoplot/{sample}NanoPlot-report.html` | NanoPlot report |

## Long reads: SNVs and indels

`{caller}` is `clair3` or `deepvariant`.

| File | Content |
| --- | --- |
| `results/tgs/variants/single_sample/{sample}/{caller}/{sample}.{caller}.vcf.gz` | Raw calls; GVCF in `{sample}.{caller}.gvcf.gz` |
| `results/tgs/variants/single_sample/{sample}/{caller}/{sample}.{caller}.norm.vcf.gz` | Normalised VCF |
| `results/tgs/variants/single_sample/{sample}/{caller}/{sample}.{caller}.stats.txt` | bcftools stats |
| `results/tgs/variants/single_sample/{sample}/deepvariant/{sample}.visual_report.html` | DeepVariant report |
| `results/tgs/variants/single_sample/{sample}/{caller}/annotated/{sample}.snpeff.vcf.gz` | snpEff annotation; table in `{sample}.snpeff.tsv` |
| `results/tgs/variants/single_sample/{sample}/{caller}/annovar/` | ANNOVAR tables and VCF |
| `results/tgs/variants/cohort/filtered/cohort.filtered.recode.vcf.gz` | Filtered cohort VCF (GLnexus or GATK) |
| `results/tgs/variants/cohort/phased/{chrom}.phased.vcf.gz` | SHAPEIT5-phased cohort VCFs |

In ANNOVAR output, `high_impact` marks exonic and splicing regions; it is not
the snpEff HIGH impact class.

## Long reads: structural variants

| File | Content |
| --- | --- |
| `results/tgs/structural_variants/single_sample/{sample}/{caller}/{sample}.{caller}.vcf` | Per-caller calls: `sniffles`, `cutesv`, `svim`, `nanovar` (ONT), `pbsv` (HiFi) |
| `results/tgs/structural_variants/single_sample/{sample}/hitsv/{sample}.hitsv.vcf.gz` | HitSV long-read calls (optional) |
| `results/tgs/structural_variants/single_sample/{sample}/merged/{sample}.merged.filtered.vcf.gz` | Per-sample merged and filtered SVs |
| `results/tgs/structural_variants/single_sample/{sample}/force_called/{sample}.force_called.vcf.gz` | cuteFC re-genotyping against the cohort candidates |
| `results/tgs/structural_variants/single_sample/{sample}/final/{sample}.sv.final.vcf.gz` | Final per-sample SVs (force-called if enabled) |
| `results/tgs/structural_variants/single_sample/{sample}/annotation/{sample}.tgs_survivor.annotsv.tsv` | AnnotSV annotation of the merged SVs (the file name is kept for compatibility, also for Truvari merging) |
| `results/tgs/structural_variants/single_sample/{sample}/annotation/{sample}.cuteFC.annotsv.tsv` | AnnotSV annotation of the force-called SVs |
| `results/tgs/structural_variants/cohort/population.merged.vcf` | Cohort SV merge, keeping only sites carried by at least one sample (`tgs_config.survivor.drop_non_carrier_sites`) |

## Hybrid and cross-platform SVs

| File | Content |
| --- | --- |
| `results/cross_platform/hitsv/{sample}/{sample}.hitsv.vcf.gz` | HitSV Hybrid calls for each pair in `samples_pairs_csv` |
| `results/annotations/hitsv/{mode}/{sample}.hitsv.annotsv.tsv` | AnnotSV annotation of HitSV calls (`hitsv.annotate`) |
| `results/cross_platform/structural_variants/{pair}/{pair}.cross_sv.merged.vcf.gz` | Union of the short- and long-read SV calls of a pair |
| `results/cross_platform/structural_variants/{pair}/{pair}.concordance.tsv` | Counts per SV type of calls found by short reads only, long reads only, or both |

In the merged VCF, `SUPP_VEC` is `11` (both), `10` (short reads only) or `01`
(long reads only). The comparison is a union merge with relaxed type matching;
it is not an accuracy measure.

| File | Content |
| --- | --- |
| `results/cross_platform/cohort_sv/cohort.cross_sv.vcf.gz` | Site-only catalogue integrating the short- and long-read cohort SV VCFs |
| `results/cross_platform/cohort_sv/cohort.cross_sv.summary.tsv` | Site counts by platform support and SV type |

Each catalogue site has `PLATFORM` (`both`, `short_only`, `long_only`), `SUPP_VEC`
(`NGS`,`TGS` order), and `NGS_AF`/`TGS_AF` with `AC`/`AN` from the cohort in which
it was found. The cohorts may contain different individuals, so the catalogue has no
genotype columns. If several sites from one cohort collapse into one catalogue
record, its AF/AC/AN come from the input site with the highest AF, not from
re-genotyping the combined site. Size/type filtering can exclude input sites;
SURVIVOR may also retain only one source ID for nearby sites in the same cohort.

## Logs and resource records

| Path | Content |
| --- | --- |
| `logs/{step}/` | Rule logs |
| `logs/benchmarks/{step}/` | Snakemake benchmark records: runtime, CPU time, memory, I/O |
| `logs/ngs_samples_*.csv`, `logs/tgs_samples_*.csv` | Samples resolved from the sample sheets and directory scans |
| `logs_cluster/` | Slurm job logs (`run_snakemake.sh` and the Slurm profiles) |

`scripts/snakemake/aggregate_benchmarks.py` summarises the benchmark records
as a table, JSON or CSV. It reads `*.txt` records only, and its
`sample_count` column counts records, not distinct samples.
