#!/bin/bash
# Script to download and prepare ANNOVAR databases for hg38
# This script helps users set up ANNOVAR annotation databases

# Exit on error
set -e

# Default settings
ANNOVAR_PATH="${1:-bin/annovar/annovar}"
DB_DIR="${2:-ref_data/annovar/humandb}"
BUILDVER="${3:-hg38}"

echo "==========================================="
echo "ANNOVAR Database Download Script"
echo "==========================================="
echo "ANNOVAR Path: $ANNOVAR_PATH"
echo "Database Directory: $DB_DIR"
echo "Build Version: $BUILDVER"
echo "==========================================="

# Check if ANNOVAR is installed
if [ ! -f "$ANNOVAR_PATH/annotate_variation.pl" ]; then
    echo "ERROR: ANNOVAR not found at $ANNOVAR_PATH"
    echo ""
    echo "Please download ANNOVAR from:"
    echo "  http://www.openbioinformatics.org/annovar/annovar_download_form.php"
    echo ""
    echo "After registration and download, extract to: $ANNOVAR_PATH"
    exit 1
fi

# Create database directory
mkdir -p "$DB_DIR"

echo ""
echo "Starting database downloads..."
echo "This may take a while depending on your internet connection."
echo ""

# Download refGene (gene-based annotation)
echo "[1/7] Downloading refGene..."
perl "$ANNOVAR_PATH/annotate_variation.pl" -buildver "$BUILDVER" \
    -downdb -webfrom annovar refGene "$DB_DIR"

# Download avsnp150 (dbSNP150)
echo "[2/7] Downloading avsnp150 (dbSNP)..."
perl "$ANNOVAR_PATH/annotate_variation.pl" -buildver "$BUILDVER" \
    -downdb -webfrom annovar avsnp150 "$DB_DIR"

# Download ClinVar
echo "[3/7] Downloading ClinVar..."
perl "$ANNOVAR_PATH/annotate_variation.pl" -buildver "$BUILDVER" \
    -downdb -webfrom annovar clinvar_20221231 "$DB_DIR"

# Download gnomAD genome
echo "[4/7] Downloading gnomAD genome (v4.1)..."
perl "$ANNOVAR_PATH/annotate_variation.pl" -buildver "$BUILDVER" \
    -downdb -webfrom annovar gnomad41_genome "$DB_DIR"

# Download dbnsfp42a (functional prediction scores)
echo "[5/7] Downloading dbNSFP (v4.2a)..."
perl "$ANNOVAR_PATH/annotate_variation.pl" -buildver "$BUILDVER" \
    -downdb -webfrom annovar dbnsfp42a "$DB_DIR"

# Download CADD scores (optional but recommended)
echo "[6/7] Downloading CADD scores..."
perl "$ANNOVAR_PATH/annotate_variation.pl" -buildver "$BUILDVER" \
    -downdb -webfrom annovar cadd13 "$DB_DIR" || echo "CADD download failed (optional)"

# Download COSMIC (optional, may require registration)
echo "[7/7] Downloading COSMIC..."
perl "$ANNOVAR_PATH/annotate_variation.pl" -buildver "$BUILDVER" \
    -downdb -webfrom annovar cosmic70 "$DB_DIR" || echo "COSMIC download failed (may require separate registration)"

echo ""
echo "==========================================="
echo "Database download completed!"
echo "==========================================="
echo ""
echo "Downloaded databases:"
ls -lh "$DB_DIR"/*.txt 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
echo ""
echo "You can now use ANNOVAR annotation in the pipeline by setting:"
echo "  tgs_variant_annotation:"
echo "    tool: \"annovar\""
echo "    annovar:"
echo "      database_dir: \"$DB_DIR\""
echo ""
