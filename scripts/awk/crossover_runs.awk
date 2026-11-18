# Ordered runs of concordant/discordant discriminating non-damage sites between two children.
# A single long run = no paternal recombination between them in the window; isolated
# single-site switches are treated as genotyping error.
# Vars: A, B = column names of the two children (e.g. A=AD143 B=AD146)
BEGIN{ FS=OFS="\t" }
NR==1{ for(i=1;i<=NF;i++) col[$i]=i; nk=NF-4; next }
$(NF)=="no" && $col[A]!="NA" && $col[B]!="NA" {
  c=""; for(i=4;i<=NF-1;i++) if($i!="NA") c=c $i
  if(c ~ /^0+$/ || c ~ /^1+$/) next
  s=($col[A]==$col[B]) ? A"="B : A"!="B
  if(s!=prev){ if(n) print start, prev, n; start=$1; n=0; prev=s }
  n++ }
END{ if(n) print start, prev, n }
