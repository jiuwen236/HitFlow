# Prepare a site-only cohort SV VCF for cross-platform matching (SURVIVOR / truvari).
# Usage: awk -v label=NGS -f prepare_sv_sites.awk sites.vcf > out.vcf
#  * ID is replaced by <label>_<n>, unique per platform; the source ID goes to INFO/SRC_ID.
#  * REF "." -> "N", <DUP:TANDEM> -> <DUP>, multi-valued SVLEN reduced to its first integer
#    (SURVIVOR 1.0.7 cannot parse these).
#  * One placeholder sample column (GT 0/1) is added; SURVIVOR needs a sample per input.
BEGIN { OFS = "\t" }
/^##/ { print; next }
/^#CHROM/ {
    print "##INFO=<ID=SRC_ID,Number=1,Type=String,Description=\"Variant ID in the source cohort VCF\">"
    print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print $1, $2, $3, $4, $5, $6, $7, $8, "FORMAT", label "_COHORT"
    next
}
{
    n++
    src = $3
    if ($4 == ".") $4 = "N"
    gsub(/<DUP:TANDEM>/, "<DUP>", $5)
    info = ""
    m = split($8, kv, ";")
    for (i = 1; i <= m; i++) {
        if (kv[i] == "" || kv[i] == ".") continue
        if (kv[i] ~ /^SVLEN=/) {
            v = substr(kv[i], 7)
            sub(/,.*/, "", v)
            if (v !~ /^-?[0-9]+$/) continue
            kv[i] = "SVLEN=" v
        }
        info = (info == "" ? kv[i] : info ";" kv[i])
    }
    info = (info == "" ? "" : info ";") "SRC_ID=" src
    print $1, $2, label "_" n, $4, $5, $6, $7, info, "GT", "0/1"
}
