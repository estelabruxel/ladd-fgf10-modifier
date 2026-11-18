#!/usr/bin/env bash
# Mean depth of the original and resequenced AD143 libraries in each target region (Table 1).
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
BED="$(dirname "$0")/../config/coverage_regions.bed"
mean_depth () { samtools depth -a -r "$1" "$2" | awk '{s+=$3; n++} END{printf "%.2f", (n? s/n : 0)}'; }

printf "region\toriginal\tresequenced\tfold_change\n" > "$OUT/AD143_coverage.tsv"
grep -v '^#' "$BED" | while IFS=$'\t' read -r C S E NAME; do
  R="${C}:$((S+1))-${E}"
  O=$(mean_depth "$R" "$AD143_BAM_ORIGINAL"); N=$(mean_depth "$R" "$AD143_BAM_RESEQ")
  awk -v n="$NAME" -v o="$O" -v r="$N" 'BEGIN{printf "%s\t%s\t%s\t%s\n", n, o, r, (o>0? sprintf("%.1f", r/o) : "NA")}'
done >> "$OUT/AD143_coverage.tsv"
column -t "$OUT/AD143_coverage.tsv"
