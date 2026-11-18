#!/usr/bin/env bash
# Joint genotyping of the 10 family members from per-sample, per-chromosome gVCFs
# (GATK HaplotypeCaller -ERC GVCF output from bcbio-nextgen; AD143 regenerated from
# the resequenced BAM). CombineGVCFs + GenotypeGVCFs per chromosome, then concatenation.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
SHEET="$(dirname "$0")/../config/gvcf_paths.tsv"
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
W="$OUT/joint"; mkdir -p "$W/combined" "$W/joint" "$W/logs"

for CHR in $CHROMS; do
  OUTV="$W/joint/family_${CHR}.joint.vcf.gz"
  [[ -s "$OUTV.tbi" ]] && { echo "$CHR: done, skipping"; continue; }
  ARGS=()
  while IFS=$'\t' read -r S P; do
    [[ -z "$S" || "$S" == \#* ]] && continue
    ARGS+=( -V "${P//\{CHR\}/$CHR}" )
  done < "$SHEET"
  $GATK --java-options "-Xmx8g" CombineGVCFs -R "$REF" -L "$CHR" "${ARGS[@]}" \
        -O "$W/combined/combined.${CHR}.g.vcf.gz" > "$W/logs/combine.${CHR}.log" 2>&1
  $GATK --java-options "-Xmx10g" GenotypeGVCFs -R "$REF" \
        -V "$W/combined/combined.${CHR}.g.vcf.gz" -O "$OUTV" > "$W/logs/genotype.${CHR}.log" 2>&1
  echo "$CHR: genotyped"
done

ls "$W"/joint/family_chr*.joint.vcf.gz | sort -V > "$W/concat.list"
bcftools concat -Oz -o "$OUT/family.joint_v2.vcf.gz" -f "$W/concat.list"
bcftools index -t "$OUT/family.joint_v2.vcf.gz"
bcftools stats "$OUT/family.joint_v2.vcf.gz" | grep -E "^SN\s+0\s+(number of samples|number of records|number of SNPs|number of indels)"
