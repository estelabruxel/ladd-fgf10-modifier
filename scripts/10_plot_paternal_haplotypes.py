#!/usr/bin/env python3
"""
Figure: paternal haplotype inheritance at the FGF10, FGFR2 and TBX4 loci.

Reads results/<LOCUS>_paternal_alleles.tsv (from 04_paternal_alleles.sh) and plots, for each
locus, one dot per discriminating site (children with a call do not all carry the same paternal
allele) at its genomic position. Colour = paternal allele relative to AD144 (haplotype "A";
AD145 is used where AD144 has no call). Hollow = C>T/G>A (damage-prone). Grey tick = no call.
Counts at the right of each row: same / other paternal allele as the reference, at
discriminating non-damage sites.

Usage: python3 10_plot_paternal_haplotypes.py [results_dir]
Output: <results_dir>/paternal_haplotypes_3loci.pdf and .png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "results")
KIDS = ["AD144", "AD145", "AD146", "AD143"]          # plotting order; AD143 last
REF_ORDER = ["AD144", "AD145"]                      # haplotype "A" reference
A_COL, B_COL, MISS = "#1a9e77", "#d95f02", "#bdbdbd"

# Annotation tracks (GRCh38). Coordinates marked VERIFY must be checked against GENCODE/VEP.
LOCI = [
    dict(file="FGF10_TAD", chrom="chr5", start=43_600_000, end=44_825_000, title="FGF10 locus",
         tracks=[("FGF10 TAD", 43_600_000, 44_825_000, "#dcd8f5"),
                 ("FGF10 deletion (maternal)", 44_300_386, 44_312_544, "#e8807a"),
                 ("FGF10 enhancer", 44_385_696, 44_386_743, "#7b6fd6")],
         marks=[(43_619_485, "rs78368856")]),
    dict(file="FGFR2_TAD", chrom="chr10", start=120_925_001, end=121_925_000, title="FGFR2 locus",
         tracks=[("FGFR2 TAD", 120_925_001, 121_925_000, "#dcd8f5"),
                 ("FGFR2 (VERIFY)", 121_478_332, 121_598_458, "#7fb8a4")],
         marks=[(121_524_165, "TTAA>T")]),
    dict(file="TBX4_region", chrom="chr17", start=61_000_000, end=61_800_000, title="TBX4 locus",
         tracks=[("TBX4 TAD", 61_348_835, 61_575_000, "#dcd8f5"),
                 ("TBX4 lung super-enhancer", 61_201_663, 61_384_701, "#7b6fd6"),
                 ("TBX4 (VERIFY)", 61_452_422, 61_485_110, "#7fb8a4")],
         marks=[]),
]


def load(path):
    rows = []
    with open(path) as f:
        head = f.readline().rstrip("\n").split("\t")
        idx = {k: head.index(k) for k in KIDS}
        for line in f:
            t = line.rstrip("\n").split("\t")
            pa = {k: (None if t[idx[k]] == "NA" else int(t[idx[k]])) for k in KIDS}
            called = {v for v in pa.values() if v is not None}
            if len(called) < 2:              # non-discriminating
                continue
            rows.append(dict(pos=int(t[0]), pa=pa, dmg=(t[-1] == "yes")))
    return rows


fig, axes = plt.subplots(len(LOCI), 1, figsize=(11, 3.4 * len(LOCI)))
for ax, L in zip(axes, LOCI):
    rows = load(os.path.join(RES, f"{L['file']}_paternal_alleles.tsv"))
    for j, kid in enumerate(KIDS):
        y = len(KIDS) - 1 - j
        ax.axhline(y, color="#f2f2f2", lw=6, zorder=0)
        same = other = 0
        for r in rows:
            ref = next((r["pa"][k] for k in REF_ORDER if r["pa"][k] is not None), None)
            v = r["pa"][kid]
            x = r["pos"] / 1e6
            if v is None or ref is None:
                ax.scatter(x, y, s=10, marker="|", color=MISS, linewidths=0.7, zorder=1)
                continue
            c = A_COL if v == ref else B_COL
            if not r["dmg"]:
                same += v == ref
                other += v != ref
            ax.scatter(x, y, s=10, marker="o", facecolors="none" if r["dmg"] else c,
                       edgecolors=c, linewidths=0.7, zorder=2)
        ax.text(1.005, y, f"{same} same / {other} other", transform=ax.get_yaxis_transform(),
                fontsize=7.5, va="center", color="#444444")
    top = len(KIDS) - 0.2
    for i, (lab, s, e, col) in enumerate(L["tracks"]):
        yb = top + 0.15 + 0.38 * i
        ax.add_patch(plt.Rectangle((s / 1e6, yb), max((e - s) / 1e6, 0.002), 0.22, color=col, lw=0))
        ax.text(L["end"] / 1e6, yb + 0.11, f"  {lab}", fontsize=7, va="center")
    for pos, lab in L["marks"]:
        ax.axvline(pos / 1e6, color="#555555", lw=0.6, ls=":", zorder=0)
        ax.text(pos / 1e6, -0.75, f" {lab}", fontsize=7, ha="left", va="top", color="#333333")
    ax.set_yticks(range(len(KIDS)))
    ax.set_yticklabels(KIDS[::-1])
    ax.get_yticklabels()[0].set_fontweight("bold")
    ax.set_ylim(-1.0, top + 0.2 + 0.38 * len(L["tracks"]))
    ax.set_xlim(L["start"] / 1e6, L["end"] / 1e6)
    ax.set_xlabel(f"{L['chrom']} position (Mb, GRCh38)", fontsize=9)
    ax.set_title(f"{L['title']}  ({len(rows)} discriminating sites)", loc="left", fontsize=10)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

handles = [plt.Line2D([], [], marker="o", ls="", color=A_COL, label="paternal haplotype A (AD144)"),
           plt.Line2D([], [], marker="o", ls="", color=B_COL, label="paternal haplotype B"),
           plt.Line2D([], [], marker="o", ls="", mfc="none", mec="#555555", label="C>T / G>A (damage-prone)"),
           plt.Line2D([], [], marker="|", ls="", color=MISS, label="no call")]
fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8, frameon=False)
fig.tight_layout(rect=(0, 0.03, 0.86, 1))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(RES, f"paternal_haplotypes_3loci.{ext}"), dpi=300)
print("Wrote", os.path.join(RES, "paternal_haplotypes_3loci.pdf/.png"))
