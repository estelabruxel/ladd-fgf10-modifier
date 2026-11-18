# Deduce the paternal allele of each child at paternally informative sites.
# Input : `bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ]\n'` (biallelic records)
# Vars  : FATHER, MOTHER, KIDS (comma-separated IDs, e.g. "AD143,AD144,AD145,AD146"), GQMIN
# Output: pos, change, mom=<GT>, paternal allele per kid (1 = ALT, 0 = REF, NA), damage (C>T/G>A)
# Informative site: father 0/1 and mother 0/0 or 1/1, both GQ >= GQMIN.
#   mother 0/0: child 0/1 -> paternal ALT; child 0/0 -> paternal REF
#   mother 1/1: child 1/1 -> paternal ALT; child 0/1 -> paternal REF
BEGIN{ FS=OFS="\t"; if(GQMIN=="") GQMIN=30; nk=split(KIDS,K,",")
       printf "pos\tchange\tmother"; for(i=1;i<=nk;i++) printf "\t%s", K[i]; print "\tdamage" }
function pat(m,ch){ if(m=="0/0"){if(ch=="0/1")return 1; if(ch=="0/0")return 0}
                    if(m=="1/1"){if(ch=="1/1")return 1; if(ch=="0/1")return 0}; return "NA" }
NR==1{ for(i=1;i<=NF;i++){ h=$i; sub(/^# ?/,"",h); sub(/^\[[0-9]+\]/,"",h)
         split(h,a,":"); s=a[1]; sub(/_.*/,"",s)
         if(a[2]=="GT") gt[s]=i; if(a[2]=="GQ") gq[s]=i }
       next }
{ d=$gt[FATHER]; m=$gt[MOTHER]; gsub(/\|/,"/",d); gsub(/\|/,"/",m)
  if(d!="0/1" || (m!="0/0" && m!="1/1")) next
  if($gq[FATHER]+0<GQMIN || $gq[MOTHER]+0<GQMIN) next
  if($4 ~ /,/ || $4=="*") next
  dmg=(length($3)==1 && length($4)==1 && (($3=="C"&&$4=="T")||($3=="G"&&$4=="A"))) ? "yes" : "no"
  line=$2 OFS $3">"$4 OFS "mom="m
  for(i=1;i<=nk;i++){ g=$gt[K[i]]; gsub(/\|/,"/",g); p[i]=pat(m,g); line=line OFS p[i] }
  print line, dmg
  # summary: the three siblings (kids 2..4) and proband (kid 1) vs kid 4
  if(p[2]!="NA"&&p[3]!="NA"&&p[4]!="NA"){ n++; if(p[2]==p[3]&&p[3]==p[4]) sib++ }
  if(p[1]!="NA"&&p[4]!="NA"&&dmg=="no"){ m14++; if(p[1]==p[4]) s14++ }
}
END{ printf "Informative sites with %s/%s/%s called: %d; all three agree: %d\n", K[2],K[3],K[4], n, sib > "/dev/stderr"
     printf "Non-damage sites with %s+%s called: %d; %s = %s: %d\n", K[1],K[4], m14, K[1],K[4], s14 > "/dev/stderr" }
