set -eo pipefail

input_vcf=$1
basename=$2
#cohort=$3
#center=$4
output_dir=$3
bcftools query  -f "[%CHROM\tCohort\t%POS\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/SR\t%SAMPLE\t%GT\n]" ${input_vcf} \
| awk 'BEGIN{OFS="\t"}{if($1~/chr[1-9]+/ && $1!~/_/) {
 svlen=$6;
 if($6<0 && $6!=".") svlen=-1*$6;
 len_bin=">=1kb"
 if(svlen<1000) len_bin="<1kb";
 if($7>0) $7="SR>=1";
 else $7="SR=0";
 print $1, $2, $3, $4, $5, len_bin, $7, $8, $9;}}' \
| sort -k1,9 \
| uniq -c \
| awk 'BEGIN{OFS="\t"}{print $2, $3, $4, $5, $6, $7, $8, $9, $1}'  > ${output_dir}/${basename}.lumpy.counts.1.txt
