# For each site, does the proband carry the ALT allele, and is it absent from every affected relative?
# Input : `bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%GQ]\n'`
# Vars  : PROBAND, AFF (comma-separated affected relatives), OUTF (per-site output)
# Status: IN_AFFECTED | ABSENT_uncertain (a relative ./. or GQ < 20) | ABSENT_all_affected
BEGIN{ FS=OFS="\t"; n=split(AFF,A,","); print "pos","change","proband_GT","status","carriers" > OUTF }
NR==1{ for(i=1;i<=NF;i++){ h=$i; sub(/^# ?/,"",h); sub(/^\[[0-9]+\]/,"",h)
         split(h,a,":"); s=a[1]; sub(/_.*/,"",s); if(a[2]=="GT")gt[s]=i; if(a[2]=="GQ")gq[s]=i }
       next }
{ g=$gt[PROBAND]; gsub(/\|/,"/",g)
  if(g!="0/1" && g!="1/1"){ notin++; next }
  car=0; unc=0; who=""
  for(i=1;i<=n;i++){ x=$gt[A[i]]; gsub(/\|/,"/",x); q=$gq[A[i]]+0
     if(x=="0/1"||x=="1/1"){ car++; who=who A[i]"," } else if(x!="0/0"||q<20) unc++ }
  st = car ? "IN_AFFECTED" : (unc ? "ABSENT_uncertain" : "ABSENT_all_affected")
  cnt[st]++
  print $2, $3">"$4, g, st, (who==""?"-":who) > OUTF }
END{ print PROBAND " does not carry ALT:", notin+0; for(s in cnt) print s":", cnt[s] }
