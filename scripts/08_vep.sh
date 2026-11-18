#!/usr/bin/env bash
# Ensembl VEP 115.1 (Docker, offline cache), one consequence per variant.
# Usage: 08_vep.sh <input.vcf> <output.tsv> [chrom]
#   If DBNSFP_DIR/dbNSFP4.1a_variant.<chrom>.bgz exists, REVEL, ClinPred, GERP++ and dbNSFP
#   phyloP fields are added (dbNSFP covers coding/splice positions only).
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
IN=$(readlink -f "$1"); OUTF=$(readlink -f "$2"); CHR=${3:-}
WD=$(dirname "$IN"); OD=$(dirname "$OUTF")

PLUGIN_ARGS=()
DBN="$DBNSFP_DIR/dbNSFP4.1a_variant.${CHR}.bgz"
if [[ -n "$CHR" && -s "$DBN" ]]; then
  PLUGIN_ARGS=( --dir_plugins /plugins
    --plugin "dbNSFP,/dbnsfp/$(basename "$DBN"),REVEL_score,phyloP100way_vertebrate,phyloP30way_mammalian,GERP++_RS,ClinPred_score,ClinPred_pred" )
fi

docker run --rm --user "$(id -u):$(id -g)" \
  -v "$VEP_CACHE":/vep_data -v "$VEP_PLUGINS":/plugins -v "$DBNSFP_DIR":/dbnsfp \
  -v "$WD":/in -v "$OD":/out \
  "$VEP_IMAGE" vep \
  --input_file "/in/$(basename "$IN")" --output_file "/out/$(basename "$OUTF")" --force_overwrite \
  --offline --cache --dir_cache /vep_data --assembly GRCh38 \
  --fasta "/vep_data/$VEP_FASTA_IN_CACHE" \
  --tab --symbol --hgvs --canonical --mane --biotype --regulatory --variant_class \
  --check_existing --af --af_gnomade --af_gnomadg --max_af \
  --pick --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length \
  --fork "$THREADS" "${PLUGIN_ARGS[@]}"

echo "VEP rows: $(grep -vc '^#' "$OUTF")"
