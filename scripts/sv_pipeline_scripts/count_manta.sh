set -eo pipefail

input_vcf=$1
basename=$2
#cohort=$3
#center=$4
output_dir=$3
bcftools query  -f "[%CHROM\tCohort\t%POS\t%FILTER\t%INFO/SVTYPE\t%SAMPLE\t%GT\n]" ${input_vcf} \
| awk 'BEGIN{OFS="\t"}{if($1~/chr[1-9]+/ && $1!~/_/ && $4=="PASS") print $0;}' \
| sort -k1,7 \
| uniq -c \
| awk 'BEGIN{OGS="\t"}{print $2, $3, $4, $5, $6, $7, $8,  $1}'  > ${output_dir}/${basename}.manta.counts.1.txt
