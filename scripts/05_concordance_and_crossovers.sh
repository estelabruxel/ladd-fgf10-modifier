#!/usr/bin/env bash
# Pairwise concordance between siblings and crossover check for each locus.
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
D="$(dirname "$0")/awk"
for F in "$OUT"/*_paternal_alleles.tsv; do
  NAME=$(basename "$F" _paternal_alleles.tsv)
  echo "=== $NAME: pairwise concordance (discriminating, non-damage) ==="
  awk -f "$D/pairwise_concordance.awk" "$F" | sort
  for B in $SIBS; do
    echo "--- runs: $PROBAND vs $B (start, state, n sites) ---"
    awk -v A="$PROBAND" -v B="$B" -f "$D/crossover_runs.awk" "$F"
  done
  echo
done | tee "$OUT/concordance_summary.txt"
