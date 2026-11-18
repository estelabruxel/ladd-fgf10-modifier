#!/usr/bin/env bash
# FGF10 TAD: variants on the paternal hap A (AD143 + unaffected sisters), absent from AD146
# and the affected aunts. Then annotation with VEP (08) and CADD/phyloP (09).
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
D="$(dirname "$0")"
REGION=chr5:43600000-44825000

bcftools view -m2 -M2 -r "$REGION" "$JOINT_VCF" -Ou |
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ\t%AD]\n' |
awk -v AUNTS="$AFFECTED_AUNTS" -f "$D/awk/hapA_filter.awk" > "$OUT/FGF10_hapA_candidates.tsv"

echo "hap A candidates: $(($(wc -l < "$OUT/FGF10_hapA_candidates.tsv")-1))"
cut -f5 "$OUT/FGF10_hapA_candidates.tsv" | tail -n +2 | sort | uniq -c

# VCF for annotation
{ echo "##fileformat=VCFv4.2"; echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  tail -n +2 "$OUT/FGF10_hapA_candidates.tsv" | awk -F'\t' -v OFS='\t' '{print $1,$2,".",$3,$4,".","PASS","."}'
} > "$OUT/FGF10_hapA_candidates.vcf"

bash "$D/08_vep.sh" "$OUT/FGF10_hapA_candidates.vcf" "$OUT/FGF10_hapA_candidates.vep.tsv" chr5
python3 "$D/09_cadd_phylop.py" --chrom chr5 \
  --input "$OUT/FGF10_hapA_candidates.tsv" --output "$OUT/FGF10_hapA_candidates.cadd_phylop.tsv"
