# Join per-site status (affected_check.awk output) with VEP tab output (gene, consequence, MAX_AF).
# VEP shifts indel positions by +1 in Location, so positions are matched within +/-1 bp.
# Usage: awk -f join_vep.awk vep.tsv status.tsv
BEGIN{ FS=OFS="\t" }
FNR==NR{ if($0 ~ /^##/) next
         if($0 ~ /^#/){ for(i=1;i<=NF;i++){ h=$i; sub(/^#/,"",h); c[h]=i }; next }
         split($c["Location"],L,"[:-]")
         ann[L[2]]=$c["SYMBOL"] OFS $c["Consequence"] OFS $c["MAX_AF"]; next }
FNR==1{ print "pos","change","status","gene","consequence","MAX_AF"; next }
$4 ~ /^ABSENT/{ x="-" OFS "-" OFS "-"
  for(d=-1; d<=1; d++) if(($1+d) in ann){ x=ann[$1+d]; break }
  print $1, $2, $4, x }
