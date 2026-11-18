#!/usr/bin/env bash
# FGFR2 TAD: variants on the paternal haplotype transmitted only to the proband.
#  1. discriminating, non-damage sites called in the proband and AD146 -> VCF
#  2. VEP annotation; count rare (gnomAD MAX_AF < 0.01)
#  3. proband carries ALT and every affected relative is 0/0 (GQ >= 20)
#  4. CADD + phyloP for the variants absent from all affected relatives
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
D="$(dirname "$0")"
PA="$OUT/FGFR2_TAD_paternal_alleles.tsv"     # from 04_paternal_alleles.sh
P="$OUT/FGFR2_TAD"

# 1. discriminating non-damage sites (columns: pos change mother AD143 AD144 AD145 AD146 damage)
{ echo "##fileformat=VCFv4.2"; echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
  awk -F'\t' -v OFS='\t' 'NR>1 && $8=="no" && $4!="NA" && $7!="NA" {
      c=""; for(i=4;i<=7;i++) if($i!="NA") c=c $i
      if(c ~ /^0+$/ || c ~ /^1+$/) next
      split($2,a,">"); print "chr10",$1,".",a[1],a[2],".","PASS","." }' "$PA"
} > "$P.discriminating.vcf"
echo "discriminating non-damage sites: $(grep -vc '^#' "$P.discriminating.vcf")"

# 2. VEP + frequency summary
bash "$D/08_vep.sh" "$P.discriminating.vcf" "$P.discriminating.vep.tsv" chr10
awk -F'\t' '/^##/{next} /^#/{for(i=1;i<=NF;i++) if($i=="MAX_AF") m=i; next}
  {seen[$1]=1; af=$m; if(af=="-"||af=="") nov[$1]=1; else if(af+0<0.01) rare[$1]=1}
  END{print "sites:",length(seen); print "rare (MAX_AF<0.01):",length(rare); print "not in gnomAD:",length(nov)}' \
  "$P.discriminating.vep.tsv"

# 3. absent from all affected relatives
grep -v '^#' "$P.discriminating.vcf" | cut -f1,2 > "$P.discriminating.pos"
bcftools view -m2 -M2 -R "$P.discriminating.pos" "$JOINT_VCF" -Ou |
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ]\n' |
awk -v PROBAND="$PROBAND" -v AFF="$AFFECTED_ALL" -v OUTF="$P.vs_affected.tsv" -f "$D/awk/affected_check.awk"

awk -f "$D/awk/join_vep.awk" "$P.discriminating.vep.tsv" "$P.vs_affected.tsv" \
  | { IFS= read -r H; echo "$H"; sort -t$'\t' -k6,6g; } > "$P.proband_only_absent_affected.tsv"
awk -F'\t' 'NR>1 && $3=="ABSENT_all_affected"' "$P.proband_only_absent_affected.tsv" | wc -l |
  xargs echo "absent from all affected:"

# 4. CADD + phyloP (confident absences only)
awk -F'\t' 'NR==1 || $3=="ABSENT_all_affected"' "$P.proband_only_absent_affected.tsv" > "$P.proband_only_confident.tsv"
python3 "$D/09_cadd_phylop.py" --chrom chr10 \
  --input "$P.proband_only_confident.tsv" --output "$P.proband_only_confident.cadd_phylop.tsv"
