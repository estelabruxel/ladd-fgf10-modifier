# Candidate hypomorphic variants tagging the paternal FGF10 haplotype carried by the proband's
# unaffected sisters (hap A) and not by her affected sibling AD146 (hap B).
# Input : `bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ\t%AD]\n'` (biallelic)
# Criteria (Methods):
#   1. father 0/1 and mother 0/0, both GQ >= 30
#   2. AD144 and AD145 0/1, GQ >= 30, depth >= 10
#   3. AD146 0/0, GQ >= 30
#   4. absent from the affected aunts (0/0 with GQ >= 20)
#   5. biallelic sites only
# AD143 is reported but not required (damaged DNA: missing or unreliable calls). Her
# heterozygous calls are flagged when the ALT fraction is outside 0.35-0.65.
# Class: STRICT = not C>T/G>A; SENSITIVITY = C>T or G>A (possible deamination artefact).
BEGIN{ FS=OFS="\t"; na=split(AUNTS,AU,",")
       print "chrom","pos","ref","alt","class","AD143_GT","AD143_AD","AD143_AF","AD143_flag",
             "AD144_GT","AD145_GT","AD146_GT" }
NR==1{ for(i=1;i<=NF;i++){ h=$i; sub(/^# ?/,"",h); sub(/^\[[0-9]+\]/,"",h)
         split(h,a,":"); s=a[1]; sub(/_.*/,"",s)
         if(a[2]=="GT") gt[s]=i; if(a[2]=="GQ") gq[s]=i; if(a[2]=="AD") ad[s]=i }
       next }
function G(s,  g){ g=$gt[s]; gsub(/\|/,"/",g); return g }
function Q(s){ return $gq[s]+0 }
function DP(s,  x,n){ n=split($ad[s],x,","); return (x[1]+0)+(n>1? x[2]+0 : 0) }
{
  if($4 ~ /,/ || $4=="*") next
  if(G("AD141")!="0/1" || G("AD142")!="0/0" || Q("AD141")<30 || Q("AD142")<30) next
  if(G("AD144")!="0/1" || G("AD145")!="0/1" || Q("AD144")<30 || Q("AD145")<30) next
  if(DP("AD144")<10 || DP("AD145")<10) next
  if(G("AD146")!="0/0" || Q("AD146")<30) next
  for(i=1;i<=na;i++) if(G(AU[i])!="0/0" || Q(AU[i])<20) next
  cls=(length($3)==1 && length($4)==1 && (($3=="C"&&$4=="T")||($3=="G"&&$4=="A"))) ? "SENSITIVITY" : "STRICT"
  g=G("AD143"); split($ad["AD143"],x,","); dp=(x[1]+0)+(x[2]+0); af=(dp? (x[2]+0)/dp : 0)
  flag = (g=="0/1" && (af<0.35 || af>0.65)) ? "AF_outside_0.35-0.65" : ((g=="0/1"||g=="1/1") ? "carrier" : "not_called_or_ref")
  print $1,$2,$3,$4,cls,g,$ad["AD143"],sprintf("%.2f",af),flag,G("AD144"),G("AD145"),G("AD146")
}
