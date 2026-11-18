# Pairwise paternal-allele concordance at discriminating, non-damage sites.
# Input: output of paternal_alleles.awk (header + pos, change, mother, kid1..kidN, damage)
# Discriminating site: the children with a call do not all carry the same paternal allele.
BEGIN{ FS=OFS="\t" }
NR==1{ nk=NF-4; for(i=1;i<=nk;i++) name[i]=$(i+3); next }
$(NF)=="no"{
  c=""; for(i=4;i<=NF-1;i++) if($i!="NA") c=c $i
  if(c ~ /^0+$/ || c ~ /^1+$/) next
  for(i=1;i<nk;i++) for(j=i+1;j<=nk;j++){
    a=$(i+3); b=$(j+3); if(a=="NA"||b=="NA") continue
    k=name[i]"-"name[j]; tot[k]++; if(a==b) eq[k]++ }
}
END{ for(k in tot) printf "%s\t%d/%d same\n", k, eq[k], tot[k] }
