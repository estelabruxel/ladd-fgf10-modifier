#!/usr/bin/env bash
# Read-based + pedigree-aware phasing (WhatsHap 2.8) of the three loci in
# config/whatshap_regions.bed, using the parents and the three siblings with high-quality DNA.
# AD143 is excluded: with her included, most variants were skipped as Mendelian conflicts
# (e.g. chr5: 1,493 skipped vs 34 without her; see README).
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
REGIONS=$(grep -v '^#' "$(dirname "$0")/../config/whatshap_regions.bed" | awk '{printf "%s%s:%d-%d", (NR>1?",":""), $1, $2+1, $3}')
W="$OUT/whatshap"; mkdir -p "$W"
SAMPLES=("$FATHER" "$MOTHER" $SIBS)
VSAMPLES=$(printf "%s${SAMPLE_SUFFIX}," "${SAMPLES[@]}"); VSAMPLES=${VSAMPLES%,}

# PED file with the VCF sample names (sex: 1 = male, 2 = female; AD146 is male)
{
  echo -e "FAM\t${FATHER}${SAMPLE_SUFFIX}\t0\t0\t1\t1"
  echo -e "FAM\t${MOTHER}${SAMPLE_SUFFIX}\t0\t0\t2\t2"
  for S in $SIBS; do
    PHENO=1; [[ "$S" == "AD146" ]] && PHENO=2
    SEX=2;   [[ "$S" == "AD146" ]] && SEX=1
    echo -e "FAM\t${S}${SAMPLE_SUFFIX}\t${FATHER}${SAMPLE_SUFFIX}\t${MOTHER}${SAMPLE_SUFFIX}\t${SEX}\t${PHENO}"
  done
} > "$W/family5.ped"

bcftools view -r "$REGIONS" -s "$VSAMPLES" -m2 -M2 -v snps,indels "$JOINT_VCF" -Oz -o "$W/family5.targets.vcf.gz"
bcftools index -t "$W/family5.targets.vcf.gz"

BAMS=(); for S in "${SAMPLES[@]}"; do BAMS+=( "${BAM_TEMPLATE//\{S\}/$S}" ); done
whatshap phase --ped "$W/family5.ped" --reference "$REF" \
  -o "$W/family5.phased.vcf.gz" "$W/family5.targets.vcf.gz" "${BAMS[@]}" 2>&1 | tee "$W/whatshap_phase.log" | tail -20
bcftools index -t "$W/family5.phased.vcf.gz"
whatshap stats "$W/family5.phased.vcf.gz" --tsv "$W/whatshap_stats.tsv"

echo "=== Proportion of heterozygous genotypes phased, per sample ==="
for S in "${SAMPLES[@]}"; do
  bcftools query -s "${S}${SAMPLE_SUFFIX}" -f '[%GT]\n' "$W/family5.phased.vcf.gz" |
  awk -v s="$S" '/^(0\|1|1\|0)$/{p++} /^(0\/1|1\/0)$/{u++} END{printf "%s\tphased=%d\tunphased=%d\t%.1f%%\n", s, p, u, (p+u? 100*p/(p+u):0)}'
done
