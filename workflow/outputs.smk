# Output file definitions and helper functions
# This file contains functions that define the expected outputs for each workflow

import os


# Helper function for joint calling outputs
def get_joint_calling_outputs():
    """Get joint calling outputs based on configuration"""
    # Only return outputs if joint calling is enabled
    if not config["workflows"]["enable_ngs_joint_calling"]:
        return []
        
    outputs = []
    
    # Always produce normalized and filtered VCF
    outputs.append(f"{JOINT_CALLING_DIR}/filtered/cohort.filtered.recode.vcf.gz")
    
    # Add per-chromosome genotyped VCFs
    outputs.extend([f"{JOINT_CALLING_DIR}/genotyped/{chromosome}.genotyped.vcf.gz"
                   for chromosome in config["chromosomes"]])
    
    # Add VQSR outputs if enabled
    if config["vqsr"]["enabled"]:
        outputs.extend([
            f"{JOINT_CALLING_DIR}/vqsr/snps.recal",
            f"{JOINT_CALLING_DIR}/vqsr/snps.tranches", 
            f"{JOINT_CALLING_DIR}/vqsr/indels.recal",
            f"{JOINT_CALLING_DIR}/vqsr/indels.tranches",
            f"{JOINT_CALLING_DIR}/vqsr/cohort.vqsr.vcf.gz"
        ])
    
    # Add phasing outputs if enabled
    if config["phasing"]["enabled"]:
        skip_chroms = config["phasing"].get("skip_chromosomes", [])
        outputs.extend([f"{JOINT_CALLING_DIR}/phased/{chromosome}.phased.vcf.gz"
                       for chromosome in config["chromosomes"]
                       if chromosome not in skip_chroms])

    return outputs


# Helper function to get basic workflow outputs
def get_basic_workflow_outputs():
    """Get outputs for workflows 1-4 (basic WGS processing)"""
    # Only return outputs if basic workflows are enabled
    if not config["workflows"]["enable_basic_ngs_workflows"]:
        return []
        
    outputs = []
    sample_ids = get_all_sample_ids()
    # Get sample IDs that have FASTQ files (not BAM-only samples)
    sample_ids_with_fastq = get_sample_ids_with_fastq()
    
    # Reference genome index and dictionary files
    outputs.extend([
        config["resolved_paths"]["reference_genome"] + ".fai",
        os.path.splitext(config["resolved_paths"]["reference_genome"])[0] + ".dict"
    ])
    
    # S03 Fastp outputs (if enabled) - cleaned FASTQ files and QC reports
    # Only for samples with FASTQ files
    if config["fastp"]["enabled"]:
        for sample_id in sample_ids_with_fastq:
            outputs.extend([
                f"{NGS_PREPROCESSING_DIR}/{sample_id}/fastp/{sample_id}.clean.R1.fq.gz",
                f"{NGS_PREPROCESSING_DIR}/{sample_id}/fastp/{sample_id}.clean.R2.fq.gz",
                f"{NGS_PREPROCESSING_DIR}/{sample_id}/fastp/{sample_id}.fastp.html",
                f"{NGS_PREPROCESSING_DIR}/{sample_id}/fastp/{sample_id}.fastp.json"
            ])
    
    # S04 FastQC outputs - only for samples with FASTQ files
    for sample_id in sample_ids_with_fastq:
        for read in ["1", "2"]:
            for ext in ["html", "zip"]:
                outputs.append(f"{NGS_PREPROCESSING_DIR}/{sample_id}/fastqc/{sample_id}_{read}_fastqc.{ext}")

    # S05 FastQC filter outputs - only for samples with FASTQ files
    outputs.extend([f"{NGS_PREPROCESSING_DIR}/{sample_id}/fastqc/{sample_id}.fastqc.report" for sample_id in sample_ids_with_fastq])
    outputs.extend([f"{NGS_PREPROCESSING_DIR}/{sample_id}/fastqc/{sample_id}.fastqc.report.json" for sample_id in sample_ids_with_fastq])

    # S07 Merge SAM outputs
    outputs.extend([f"{NGS_ALIGNMENT_DIR}/{sample_id}/bam/{sample_id}.dup.sorted.bam" for sample_id in sample_ids])
    outputs.extend([f"{NGS_ALIGNMENT_DIR}/{sample_id}/bam/{sample_id}.dup.sorted.bai" for sample_id in sample_ids])
    # Split and discordant BAMs are only generated during alignment from FASTQ
    outputs.extend([f"{NGS_ALIGNMENT_DIR}/{sample_id}/bam/{sample_id}.{type}.bam" for sample_id in sample_ids_with_fastq for type in ["split", "discordant"]])
    outputs.extend([f"{NGS_ALIGNMENT_DIR}/{sample_id}/bam/{sample_id}.{type}.bai" for sample_id in sample_ids_with_fastq for type in ["split", "discordant"]])

    # S11 GATK BQSR outputs
    outputs.extend([f"{NGS_ALIGNMENT_DIR}/{sample_id}/bam/{sample_id}.dup.sorted.recal.bam" for sample_id in sample_ids])
    outputs.extend([f"{NGS_ALIGNMENT_DIR}/{sample_id}/bam/{sample_id}.dup.sorted.recal.bai" for sample_id in sample_ids])
              
    # S14 Qualimap BamQC outputs
    outputs.extend([f"{NGS_ALIGNMENT_DIR}/{sample_id}/qc/qualimap/qualimapReport.html" for sample_id in sample_ids])
    
    # S15 Qualimap Filter outputs
    outputs.extend([f"{NGS_ALIGNMENT_DIR}/{sample_id}/qc/{sample_id}.qualimap.report" for sample_id in sample_ids])
    
    # S16 & S17 Verify BAM ID outputs (only if enabled)
    if config.get("verifybamid2", {}).get("enabled", True):
        outputs.extend([f"{NGS_ALIGNMENT_DIR}/{sample_id}/qc/verifybamid/{sample_id}.selfSM" for sample_id in sample_ids])
        outputs.extend([f"{NGS_ALIGNMENT_DIR}/{sample_id}/qc/{sample_id}.verifybamid.report" for sample_id in sample_ids])
    
    # S19 GATK HaplotypeCaller outputs
    outputs.extend([f"{NGS_VARIANTS_SINGLE_DIR}/{sample_id}/gvcf/{sample_id}.{chrom}.hc.g.vcf.bgz" 
                   for sample_id in sample_ids for chrom in config["chromosomes"]])
    
    # Single-sample variant filtering and annotation outputs
    if config["workflows"].get("enable_ngs_single_sample_filtering", False):
        # Filtered and normalized VCF
        outputs.extend([f"{NGS_VARIANTS_SINGLE_DIR}/{sample_id}/variant_calling/normalized/{sample_id}.filtered.normalized.vcf.gz"
                       for sample_id in sample_ids])
        
        # Annotation outputs (based on selected tool)
        annotation_tool = config.get("ngs_single_sample_variants", {}).get("annotation", {}).get("tool", "funcotator")
        if annotation_tool == "funcotator":
            outputs.extend([f"{NGS_VARIANTS_SINGLE_DIR}/{sample_id}/variant_calling/annotated/{sample_id}.funcotator.vcf.gz"
                           for sample_id in sample_ids])
            outputs.extend([f"{NGS_VARIANTS_SINGLE_DIR}/{sample_id}/variant_calling/annotated/{sample_id}.funcotator.tsv"
                           for sample_id in sample_ids])
        elif annotation_tool == "snpeff":
            outputs.extend([f"{NGS_VARIANTS_SINGLE_DIR}/{sample_id}/variant_calling/annotated/{sample_id}.snpeff.vcf.gz"
                           for sample_id in sample_ids])
            outputs.extend([f"{NGS_VARIANTS_SINGLE_DIR}/{sample_id}/variant_calling/annotated/{sample_id}.snpeff.tsv"
                           for sample_id in sample_ids])
            outputs.extend([f"{NGS_VARIANTS_SINGLE_DIR}/{sample_id}/variant_calling/annotated/{sample_id}.snpeff_summary.html"
                           for sample_id in sample_ids])
    
    return outputs


# Helper function for SV detection outputs
def get_sv_detection_outputs():
    """Get SV detection outputs based on configuration"""
    # Only return outputs if SV detection is enabled
    if not config["workflows"]["enable_tgs_sv_detection"]:
        return []
        
    outputs = []
    ont_samples = get_ont_samples()
    
    if not ont_samples:  # No ONT samples found
        return []
    
    sv_dir = TGS_SV_SINGLE_DIR
    sv_cohort_dir = TGS_SV_COHORT_DIR
    
    # Get SV detection switches
    _sv_cfg = config["tgs_config"].get("sv_detection", {})
    enable_merge = _sv_cfg.get("merge_per_sample", True)
    enable_force = _sv_cfg.get("force_calling", True)
    enable_pop = _sv_cfg.get("population_merge", True)
    
    sv_filter_enabled = config["tgs_config"].get("sv_filter", {}).get("enabled", True)

    # A caller-only run has no implicit QC, merge, annotation or final-merge targets.
    if not (enable_merge or enable_force or enable_pop):
        from types import SimpleNamespace
        for sample in ont_samples:
            outputs.extend(_merged_vcf_inputs(SimpleNamespace(sample=sample)))
        return outputs

    # Individual sample outputs (based on switches)
    for sample in ont_samples:
        # QC output resides in preprocessing directory
        outputs.append(f"{TGS_PREPROCESSING_DIR}/{sample}/qc/nanoplot/{sample}NanoPlot-report.html")
        
        # Per-sample merged VCF (if merge is enabled)
        if enable_merge:
            outputs.append(f"{sv_dir}/{sample}/merged/{sample}.merged.vcf")
            if sv_filter_enabled:
                outputs.append(f"{sv_dir}/{sample}/merged/{sample}.merged.filtered.vcf.gz")
                outputs.append(f"{sv_dir}/{sample}/merged/{sample}.merged.filtered.vcf.gz.tbi")
            if config.get("annotsv", {}).get("enabled", True):
                outputs.append(f"{sv_dir}/{sample}/annotation/{sample}.tgs_survivor.annotsv.tsv")
        
        # Force-called VCF (if force calling is enabled)
        if enable_force:
            outputs.append(f"{sv_dir}/{sample}/force_called/{sample}.force_called.vcf")
            outputs.append(f"{sv_dir}/{sample}/force_called/{sample}.force_called.vcf.gz")
            outputs.append(f"{sv_dir}/{sample}/force_called/{sample}.force_called.vcf.gz.tbi")
            if config.get("annotsv", {}).get("enabled", True):
                outputs.append(f"{sv_dir}/{sample}/annotation/{sample}.cuteFC.annotsv.tsv")
        
        # Final deliverable (force-called if available, otherwise filtered merge)
        outputs.append(f"{sv_dir}/{sample}/final/{sample}.sv.final.vcf.gz")
        outputs.append(f"{sv_dir}/{sample}/final/{sample}.sv.final.vcf.gz.tbi")
    
    # Population-level output (if population merge is enabled)
    if enable_pop:
        outputs.append(f"{sv_cohort_dir}/population.merged.vcf")
    
    return outputs


def get_tgs_variant_calling_outputs_wrapper():
    """Wrapper function for TGS variant calling outputs (Clair3 or DeepVariant)"""
    if not config["workflows"].get("enable_tgs_variant_calling", False):
        return []
    
    # Use helper function that supports "Default" auto-selection based on platform
    selected_tool = get_tgs_variant_calling_tool()
    
    try:
        if selected_tool == "clair3":
            return get_tgs_variant_calling_outputs()
        elif selected_tool == "deepvariant":
            return get_tgs_deepvariant_outputs()
        else:
            # Default to Clair3 if tool is not recognized
            return get_tgs_variant_calling_outputs()
    except NameError:
        return []


def get_ngs_sv_outputs():
    """Get NGS SV detection outputs"""
    if not config["workflows"].get("enable_ngs_sv_detection", False):
        return []
    try:
        return ngs_sv_final_outputs()
    except NameError:
        return []


def get_tgs_annotation_outputs_wrapper():
    """Wrapper to get TGS annotation outputs based on selected annotation tool."""
    if not config["workflows"].get("enable_tgs_variant_annotation", False):
        return []

    # Get selected annotation tool
    tgs_annotation_tool = config.get("tgs_variant_annotation", {}).get("tool", "snpeff").strip().lower()

    try:
        if tgs_annotation_tool == "snpeff":
            return get_tgs_annotation_outputs()  # From tgs_snpeff.smk
        elif tgs_annotation_tool == "annovar":
            return get_tgs_annovar_outputs()  # From tgs_annovar.smk
        else:
            # Fallback to snpeff if unrecognized tool
            return get_tgs_annotation_outputs()
    except NameError:
        return []


def get_cross_platform_sv_outputs():
    """Get cross-platform (NGS+TGS) SV integration outputs.

    For each validated sample pair, produce the merged SV VCF (carrying SUPP_VEC)
    and the per-pair concordance summary. Returns [] when the feature is disabled
    or no valid pairs are available.
    """
    if not config["workflows"].get("enable_cross_platform_sv", False):
        return []

    try:
        pair_ids = get_cross_platform_pair_ids()
    except NameError:
        return []

    if not pair_ids:
        return []

    cross_dir = (
        config["output"].get("cross_platform", {}).get("sv")
        if isinstance(config["output"].get("cross_platform"), dict)
        else None
    ) or os.path.join(config["output"]["root"], "cross_platform", "structural_variants")

    outputs = []
    for pair_id in pair_ids:
        outputs.extend([
            f"{cross_dir}/{pair_id}/{pair_id}.cross_sv.merged.vcf.gz",
            f"{cross_dir}/{pair_id}/{pair_id}.cross_sv.merged.vcf.gz.tbi",
            f"{cross_dir}/{pair_id}/{pair_id}.concordance.tsv",
        ])
    return outputs


def get_hitsv_outputs():
    """Expose independently enabled HitSV callsets, including Hybrid pairs."""
    return hitsv_outputs() if "hitsv_outputs" in globals() else []
