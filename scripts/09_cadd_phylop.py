#!/usr/bin/env python3
"""
Add CADD v1.7 PHRED and UCSC phyloP (100-vertebrate, 30-mammal) scores to a variant table.

Input : tab-separated file with a header containing 'pos' and either 'change' (REF>ALT)
        or 'ref' and 'alt'.
Output: the same table with CADD_PHRED, phyloP100way, phyloP30way appended.

CADD: precomputed GRCh38 files (SNVs; gnomAD v4.0 indels), queried with tabix; chromosomes
      have no 'chr' prefix. Indels absent from the indel file get NA.
phyloP: SNVs = score at the base; indels = maximum over the reference bases spanned.
Paths come from --cadd-snv/--cadd-indel/--phylop100/--phylop30 or the environment
(CADD_SNV, CADD_INDEL, PHYLOP100, PHYLOP30, exported by config/config.sh).
"""
import argparse
import csv
import os
import subprocess

import pyBigWig

ap = argparse.ArgumentParser()
ap.add_argument("--chrom", required=True, help="e.g. chr10")
ap.add_argument("--input", required=True)
ap.add_argument("--output", required=True)
ap.add_argument("--cadd-snv", default=os.environ.get("CADD_SNV"))
ap.add_argument("--cadd-indel", default=os.environ.get("CADD_INDEL"))
ap.add_argument("--phylop100", default=os.environ.get("PHYLOP100"))
ap.add_argument("--phylop30", default=os.environ.get("PHYLOP30"))
a = ap.parse_args()

bw100 = pyBigWig.open(a.phylop100)
bw30 = pyBigWig.open(a.phylop30)
cadd_chrom = a.chrom.replace("chr", "")


def cadd(pos, ref, alt):
    f = a.cadd_snv if len(ref) == 1 and len(alt) == 1 else a.cadd_indel
    out = subprocess.run(["tabix", f, f"{cadd_chrom}:{pos}-{pos}"],
                         capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        t = line.split("\t")  # chrom pos ref alt raw phred
        if int(t[1]) == pos and t[2] == ref and t[3] == alt:
            return t[5]
    return "NA"


def phylo(bw, pos, ref):
    vals = [v for v in bw.values(a.chrom, pos - 1, pos - 1 + len(ref)) if v == v]
    return f"{max(vals):.3f}" if vals else "NA"


n = 0
with open(a.input) as fin, open(a.output, "w", newline="") as fout:
    rd = csv.DictReader(fin, delimiter="\t")
    cols = rd.fieldnames + ["CADD_PHRED", "phyloP100way", "phyloP30way"]
    wr = csv.DictWriter(fout, fieldnames=cols, delimiter="\t")
    wr.writeheader()
    for row in rd:
        pos = int(row["pos"])
        ref, alt = row["change"].split(">") if "change" in row else (row["ref"], row["alt"])
        row["CADD_PHRED"] = cadd(pos, ref, alt)
        row["phyloP100way"] = phylo(bw100, pos, ref)
        row["phyloP30way"] = phylo(bw30, pos, ref)
        wr.writerow(row)
        n += 1
print(f"annotated: {n} -> {a.output}")
