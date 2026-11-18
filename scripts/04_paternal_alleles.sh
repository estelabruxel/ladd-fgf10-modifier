#!/usr/bin/env bash
# Paternal allele of each child at informative sites in the three loci (joint VCF, biallelic).
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
AWK="$(dirname "$0")/awk/paternal_alleles.awk"
KIDS="${PROBAND},$(echo $SIBS | tr ' ' ',')"

grep -v '^#' "$(dirname "$0")/../config/loci.bed" | while IFS=$'\t' read -r C S E NAME; do
  echo "=== $NAME ($C:$((S+1))-$E) ==="
  bcftools view -m2 -M2 -r "$C:$((S+1))-$E" "$JOINT_VCF" -Ou |
  bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ]\n' |
  awk -v FATHER="$FATHER" -v MOTHER="$MOTHER" -v KIDS="$KIDS" -v GQMIN=30 -f "$AWK" \
    > "$OUT/${NAME}_paternal_alleles.tsv"
done
