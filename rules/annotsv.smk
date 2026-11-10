# SV Annotation rules using AnnotSV

# Per-rule conda environment
CONDA_ANNOTSV_ENV = "../envs/sv_annotation/annotsv.yaml"
CONTAINER_ANNOTSV = config.get("containers", {}).get("annotsv", "")

# Directories
# Note: This file works with both NGS and TGS, using specific directory variables
LOGS_DIR = config["output"]["logs"]

# Configs for AnnotSV
_ANNOTSV_CFG = config.get("annotsv", {})
ANNOTSV_GENOME_BUILD = _ANNOTSV_CFG.get("genome_build", "GRCh38")
ANNOTSV_MIN_SIZE = int(_ANNOTSV_CFG.get("sv_min_size", 30))
NORMALIZE_INV_SCRIPT = workflow_path("scripts/sv_pipeline_scripts/normalize_inv_for_annotsv.py")
ANNOTSV_ANNOT_DIR = _ANNOTSV_CFG.get("annotations_dir", "ref_data/annotsv_annotations")  # directory where AnnotSV annotation data will (be) downloaded

# --------------- NGS: annotate per-sample merged SVs ---------------

rule ngs_annotsv_per_sample:
    input:
        merged_vcf = f"{NGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.ngs_sv.merged.vcf.gz"
    output:
        tsv = f"{NGS_SV_SINGLE_DIR}/{{sample}}/annotation/{{sample}}.ngs.annotsv.tsv"
    threads: 4
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    conda:
        CONDA_ANNOTSV_ENV
    container:
        CONTAINER_ANNOTSV
    benchmark:
        f"{LOGS_DIR}/benchmarks/annotsv/ngs/{{sample}}.tsv"
    log:
        f"{LOGS_DIR}/annotsv/ngs/{{sample}}.log"
    shell:
        """
        # AnnotSV treats helper-script stderr as fatal, so use a locale provided by the image.
        export LANG=C LC_ALL=C
        mkdir -p $(dirname {output.tsv}) {ANNOTSV_ANNOT_DIR} "${{TMPDIR:-/tmp}}"
        ANNOTSV_WORK_DIR=$(mktemp -d "${{TMPDIR:-/tmp}}/annotsv.{wildcards.sample}.ngs.XXXXXX")
        ANNOTSV_ANNOT_DIR_ABS=$(readlink -f "{ANNOTSV_ANNOT_DIR}")
        trap 'rm -rf "$ANNOTSV_WORK_DIR"' EXIT
        # Prepare a clean VCF for AnnotSV (keep headers, filter malformed data records)
        TMP_VCF_RAW="$ANNOTSV_WORK_DIR/input.raw.vcf"
        TMP_VCF="$ANNOTSV_WORK_DIR/input.vcf"
        bcftools view -Ov {input.merged_vcf} > "$TMP_VCF_RAW" 2>> {log}
        # Keep all header lines (## and #CHROM), but filter data records: must have valid CHROM, numeric POS, valid REF/ALT
        awk 'BEGIN{{OFS="	"}} /^##/ {{print; next}} /^#CHROM/ {{print; next}} {{if (NF >= 5 && $1 !~ /^$/ && $2 ~ /^[0-9]+$/ && $4 !~ /^$/ && $5 !~ /^$/) print}}' "$TMP_VCF_RAW" | sed 's/\r$//' > "$TMP_VCF" || true
        # Check if the filtered VCF has any data records (not just headers) - use || true for strict mode, trim whitespace from wc output
        DATA_LINES=$(grep -v "^#" "$TMP_VCF" 2>/dev/null | wc -l | tr -d '[:space:]' || echo "0")
        # Normalize sequence-resolved inversions (REF/revcomp-ALT, no END/SVLEN) into symbolic
        # <INV> with END/SVLEN so AnnotSV can size them; otherwise it derives length 0 and aborts.
        python3 {NORMALIZE_INV_SCRIPT} --input "$TMP_VCF" --output "$TMP_VCF.norm" 2>> {log}
        mv "$TMP_VCF.norm" "$TMP_VCF"
        if [ "$DATA_LINES" -eq 0 ]; then
            echo "Warning: No valid SV records found in input VCF after filtering. Creating empty output." >> {log}
            touch {output.tsv}
        else
            # AnnotSV appends .tsv; its .unannotated.tsv file is only a rejection report.
            BASE_NAME=$(basename {output.tsv} .tsv)
            AnnotSV \
                -SVinputFile "$TMP_VCF" \
                -outputDir "$ANNOTSV_WORK_DIR" \
                -outputFile "$BASE_NAME" \
                -annotationsDir "$ANNOTSV_ANNOT_DIR_ABS" \
                -genomeBuild {ANNOTSV_GENOME_BUILD} \
                -SVminSize {ANNOTSV_MIN_SIZE} >> {log} 2>&1
            
            ANNOTSV_RESULT="$ANNOTSV_WORK_DIR/${{BASE_NAME}}.tsv"
            if [ ! -s "$ANNOTSV_RESULT" ]; then
                echo "Error: AnnotSV completed without a non-empty TSV output" >> {log}
                exit 1
            fi
            cp "$ANNOTSV_RESULT" {output.tsv}
        fi
        """

# --------------- TGS: annotate per-sample SURVIVOR merged SVs ---------------

rule tgs_annotsv_after_survivor:
    input:
        merged_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/merged/{{sample}}.merged.filtered.vcf.gz"
    output:
        tsv = f"{TGS_SV_SINGLE_DIR}/{{sample}}/annotation/{{sample}}.tgs_survivor.annotsv.tsv"
    threads: 4
    resources:
        mem_mb = 8192,
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    conda:
        CONDA_ANNOTSV_ENV
    container:
        CONTAINER_ANNOTSV
    benchmark:
        f"{LOGS_DIR}/benchmarks/annotsv/tgs_survivor/{{sample}}.tsv"
    log:
        f"{LOGS_DIR}/annotsv/tgs_survivor/{{sample}}.log"
    shell:
        """
        # AnnotSV treats helper-script stderr as fatal, so use a locale provided by the image.
        export LANG=C LC_ALL=C
        mkdir -p $(dirname {output.tsv}) {ANNOTSV_ANNOT_DIR} "${{TMPDIR:-/tmp}}"
        ANNOTSV_WORK_DIR=$(mktemp -d "${{TMPDIR:-/tmp}}/annotsv.{wildcards.sample}.tgs-survivor.XXXXXX")
        ANNOTSV_ANNOT_DIR_ABS=$(readlink -f "{ANNOTSV_ANNOT_DIR}")
        trap 'rm -rf "$ANNOTSV_WORK_DIR"' EXIT
        # Ensure uncompressed and clean VCF for AnnotSV
        TMP_VCF_RAW="$ANNOTSV_WORK_DIR/input.raw.vcf"
        TMP_VCF="$ANNOTSV_WORK_DIR/input.vcf"
        bcftools view -Ov {input.merged_vcf} > "$TMP_VCF_RAW" 2>> {log}
        # Keep all header lines (## and #CHROM), but filter data records: must have valid CHROM, numeric POS, valid REF/ALT
        awk 'BEGIN{{OFS="	"}} /^##/ {{print; next}} /^#CHROM/ {{print; next}} {{if (NF >= 5 && $1 !~ /^$/ && $2 ~ /^[0-9]+$/ && $4 !~ /^$/ && $5 !~ /^$/) print}}' "$TMP_VCF_RAW" | sed 's/\r$//' > "$TMP_VCF" || true
        # Check if the filtered VCF has any data records (not just headers) - use || true for strict mode, trim whitespace from wc output
        DATA_LINES=$(grep -v "^#" "$TMP_VCF" 2>/dev/null | wc -l | tr -d '[:space:]' || echo "0")
        # Normalize sequence-resolved inversions (REF/revcomp-ALT, no END/SVLEN) into symbolic
        # <INV> with END/SVLEN so AnnotSV can size them; otherwise it derives length 0 and aborts.
        python3 {NORMALIZE_INV_SCRIPT} --input "$TMP_VCF" --output "$TMP_VCF.norm" 2>> {log}
        mv "$TMP_VCF.norm" "$TMP_VCF"
        if [ "$DATA_LINES" -eq 0 ]; then
            echo "Warning: No valid SV records found in input VCF after filtering. Creating empty output." >> {log}
            touch {output.tsv}
        else
            # AnnotSV appends .tsv; its .unannotated.tsv file is only a rejection report.
            BASE_NAME=$(basename {output.tsv} .tsv)
            AnnotSV \
                -SVinputFile "$TMP_VCF" \
                -outputDir "$ANNOTSV_WORK_DIR" \
                -outputFile "$BASE_NAME" \
                -annotationsDir "$ANNOTSV_ANNOT_DIR_ABS" \
                -genomeBuild {ANNOTSV_GENOME_BUILD} \
                -SVminSize {ANNOTSV_MIN_SIZE} >> {log} 2>&1
            
            ANNOTSV_RESULT="$ANNOTSV_WORK_DIR/${{BASE_NAME}}.tsv"
            if [ ! -s "$ANNOTSV_RESULT" ]; then
                echo "Error: AnnotSV completed without a non-empty TSV output" >> {log}
                exit 1
            fi
            cp "$ANNOTSV_RESULT" {output.tsv}
        fi
        """

# --------------- TGS: annotate per-sample cuteFC force-called SVs ---------------

rule tgs_annotsv_after_cuteFC:
    input:
        force_vcf = f"{TGS_SV_SINGLE_DIR}/{{sample}}/force_called/{{sample}}.force_called.vcf"
    output:
        tsv = f"{TGS_SV_SINGLE_DIR}/{{sample}}/annotation/{{sample}}.cuteFC.annotsv.tsv"
    threads: 4
    resources:
        mem_mb = int(_ANNOTSV_CFG.get("cutefc_mem_mb", 16384)),
        partition = config["cluster"]["partition"],
        time = config["cluster"]["time"]
    conda:
        CONDA_ANNOTSV_ENV
    container:
        CONTAINER_ANNOTSV
    benchmark:
        f"{LOGS_DIR}/benchmarks/annotsv/tgs_cuteFC/{{sample}}.tsv"
    log:
        f"{LOGS_DIR}/annotsv/tgs_cuteFC/{{sample}}.log"
    shell:
        """
        # AnnotSV treats helper-script stderr as fatal, so use a locale provided by the image.
        export LANG=C LC_ALL=C
        mkdir -p $(dirname {output.tsv}) {ANNOTSV_ANNOT_DIR} "${{TMPDIR:-/tmp}}"
        ANNOTSV_WORK_DIR=$(mktemp -d "${{TMPDIR:-/tmp}}/annotsv.{wildcards.sample}.cutefc.XXXXXX")
        ANNOTSV_ANNOT_DIR_ABS=$(readlink -f "{ANNOTSV_ANNOT_DIR}")
        trap 'rm -rf "$ANNOTSV_WORK_DIR"' EXIT
        # Ensure uncompressed and clean VCF for AnnotSV
        TMP_VCF_RAW="$ANNOTSV_WORK_DIR/input.raw.vcf"
        TMP_VCF="$ANNOTSV_WORK_DIR/input.vcf"
        bcftools view -Ov {input.force_vcf} > "$TMP_VCF_RAW" 2>> {log}
        # Keep all header lines (## and #CHROM), but filter data records: must have valid CHROM,
        # numeric POS, valid REF/ALT, AND drop translocations/breakends (SVTYPE=TRA/BND). The
        # force-called set carries ~175k TRA records (breakend pairs) that AnnotSV cannot size and
        # that blow up its runtime; they are not annotated on the merged branch either.
        awk 'BEGIN{{OFS="	"}} /^##/ {{print; next}} /^#CHROM/ {{print; next}} {{if (NF >= 5 && $1 !~ /^$/ && $2 ~ /^[0-9]+$/ && $4 !~ /^$/ && $5 !~ /^$/ && $8 !~ /SVTYPE=TRA/ && $8 !~ /SVTYPE=BND/) print}}' "$TMP_VCF_RAW" | sed 's/\r$//' > "$TMP_VCF" || true
        # Check if the filtered VCF has any data records (not just headers) - use || true for strict mode, trim whitespace from wc output
        DATA_LINES=$(grep -v "^#" "$TMP_VCF" 2>/dev/null | wc -l | tr -d '[:space:]' || echo "0")
        # Normalize sequence-resolved inversions (REF/revcomp-ALT, no END/SVLEN) into symbolic
        # <INV> with END/SVLEN so AnnotSV can size them; otherwise it derives length 0 and aborts.
        python3 {NORMALIZE_INV_SCRIPT} --input "$TMP_VCF" --output "$TMP_VCF.norm" 2>> {log}
        mv "$TMP_VCF.norm" "$TMP_VCF"
        if [ "$DATA_LINES" -eq 0 ]; then
            echo "Warning: No valid SV records found in input VCF after filtering. Creating empty output." >> {log}
            touch {output.tsv}
        else
            # AnnotSV appends .tsv; its .unannotated.tsv file is only a rejection report.
            BASE_NAME=$(basename {output.tsv} .tsv)
            AnnotSV \
                -SVinputFile "$TMP_VCF" \
                -outputDir "$ANNOTSV_WORK_DIR" \
                -outputFile "$BASE_NAME" \
                -annotationsDir "$ANNOTSV_ANNOT_DIR_ABS" \
                -genomeBuild {ANNOTSV_GENOME_BUILD} \
                -SVminSize {ANNOTSV_MIN_SIZE} >> {log} 2>&1
            
            ANNOTSV_RESULT="$ANNOTSV_WORK_DIR/${{BASE_NAME}}.tsv"
            if [ ! -s "$ANNOTSV_RESULT" ]; then
                echo "Error: AnnotSV completed without a non-empty TSV output" >> {log}
                exit 1
            fi
            cp "$ANNOTSV_RESULT" {output.tsv}
        fi
        """ 
