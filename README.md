# Paternal haplotype analysis of a non-coding modifier in trans to an *FGF10* deletion

Code for the whole-genome sequencing (WGS) analyses in:

> [Authors]. *[Title]*. [Journal, year]. doi: [DOI]

A 12 kb heterozygous deletion of the last two exons of *FGF10* segregates in a family with
lacrimo-auriculo-dento-digital (LADD) syndrome (Wade et al. 2022). All deletion carriers survived
except a neonate (AD143) who died of acinar dysplasia. This repository contains the pipeline used
to (1) reconstruct which paternal haplotype each sibling inherited at the *FGF10*, *FGFR2* and *TBX4*
loci, and (2) select and annotate candidate non-coding variants on the paternal chromosome in trans
to the maternal deletion.

## Family and sample IDs

| ID | Relationship | Phenotype | Notes |
|---|---|---|---|
| AD141 | Father | Unaffected | |
| AD142 | Mother | LADD, *FGF10* deletion | |
| AD143 | Proband | Lethal acinar dysplasia, *FGF10* deletion | Degraded DNA (C>T / G>A deamination); resequenced |
| AD144, AD145 | Sisters | Unaffected | |
| AD146 | Brother | LADD, *FGF10* deletion | |
| AD137, AD138, AD139, AD140 | Maternal relatives | LADD, *FGF10* deletion | "affected aunts" in the scripts |

## Data availability

Participant sequencing data (BAM, gVCF, VCF) are **not** included in this repository.
`.gitignore` blocks genomic data files from being committed. [Describe controlled-access deposition / request procedure.]

## Requirements

| Tool | Version used |
|---|---|
| GATK (HaplotypeCaller, CombineGVCFs, GenotypeGVCFs) | v4.6.2.0 |
| bcftools / samtools / htslib (tabix) | [1.23.1] |
| WhatsHap | 2.8 |
| Ensembl VEP (Docker `ensemblorg/ensembl-vep:release_115.1`) with the GRCh38 115 cache | 115.1 |
| CADD precomputed scores (GRCh38): `whole_genome_SNVs.tsv.gz`, `gnomad.genomes.r4.0.indel.tsv.gz` | v1.7 |
| UCSC phyloP bigWigs: `hg38.phyloP100way.bw`, `hg38.phyloP30way.bw` | hg38 |
| dbNSFP (optional VEP plugin) | 4.1a |
| Python 3 with `pyBigWig`, `matplotlib` | |

`scripts/00_versions.sh` prints the installed versions.

## Setup

1. Edit `config/config.sh` (reference, joint VCF, BAMs, VEP cache, CADD and phyloP paths, sample suffix).
2. Edit `config/gvcf_paths.tsv` (only needed for step 01).
3. Check `config/loci.bed` and `config/coverage_regions.bed`. Entries marked `CHECK` need verification.

Sample names in the VCF are expected to be `<ID><SAMPLE_SUFFIX>` (e.g. `AD141_sort_bam`); the awk
scripts strip everything after the first underscore.

## Workflow

| Step | Script | What it does | Main output (`results/`) |
|---|---|---|---|
| 01 | `01_joint_genotyping.sh` | CombineGVCFs + GenotypeGVCFs per chromosome for the 10 samples; `bcftools concat` | `family.joint_v3.vcf.gz` |
| 02 | `02_coverage.sh` | Mean depth of the original vs resequenced AD143 libraries in target regions | `AD143_coverage.tsv` |
| 03 | `03_whatshap_phasing.sh` | WhatsHap read-based + pedigree phasing of the three loci (`config/whatshap_regions.bed`) with the parents, AD144, AD145, AD146 (AD143 excluded) | `whatshap/` |
| 04 | `04_paternal_alleles.sh` | Paternal allele of each child at informative sites in the three loci | `<LOCUS>_paternal_alleles.tsv` |
| 05 | `05_concordance_and_crossovers.sh` | Pairwise sibling concordance and crossover runs | `concordance_summary.txt` |
| 06 | `06_fgf10_hapA_candidates.sh` | *FGF10* TAD variants tagging hap A (proband + unaffected sisters), absent from AD146 and affected aunts; VEP, CADD, phyloP | `FGF10_hapA_candidates.*` |
| 07 | `07_fgfr2_proband_only.sh` | *FGFR2* TAD variants on the paternal haplotype transmitted only to AD143, absent from all affected relatives; VEP, CADD, phyloP | `FGFR2_TAD.*` |
| 08 | `08_vep.sh` | VEP wrapper used by 06 and 07 | |
| 09 | `09_cadd_phylop.py` | CADD PHRED + phyloP100way/30way for a variant table | |
| 10 | `10_plot_paternal_haplotypes.py` | Three-panel figure of paternal haplotypes at the three loci | `paternal_haplotypes_3loci.pdf/.png` |

Run in order, e.g.:

```bash
bash scripts/00_versions.sh
bash scripts/04_paternal_alleles.sh
bash scripts/05_concordance_and_crossovers.sh
bash scripts/06_fgf10_hapA_candidates.sh
bash scripts/07_fgfr2_proband_only.sh
python3 scripts/10_plot_paternal_haplotypes.py
```

## Method details

**Informative site.** Father 0/1 and mother homozygous (0/0 or 1/1), both GQ ≥ 30; biallelic records only.
The paternal allele of each child is then deduced from her genotype:

| Mother | Child | Paternal allele |
|---|---|---|
| 0/0 | 0/1 | ALT |
| 0/0 | 0/0 | REF |
| 1/1 | 1/1 | ALT |
| 1/1 | 0/1 | REF |

**Discriminating site.** The children with a genotype call do not all carry the same paternal allele.
Sites where all children agree cannot distinguish the two paternal chromosomes (and usually reflect a
paternal genotyping error) and are excluded from concordance.

**Damage-aware filtering.** AD143's DNA shows cytosine deamination, so C>T and G>A substitutions are
excluded from concordance and crossover analyses and are reported separately (sensitivity set) in
candidate lists. AD143 is excluded from WhatsHap phasing; her paternal alleles are deduced directly
from her genotypes and compared with those of her siblings.

**Why AD143 is excluded from WhatsHap.** Including her made WhatsHap discard most variants as
Mendelian conflicts (WhatsHap 2.8 logs):

| Locus | 5 samples (AD143 excluded): usable / skipped | 6 samples, original AD143 BAM (`whatshap_phase.log`) | 6 samples, resequenced AD143 BAM (`whatshap6_v3.log`) |
|---|---|---|---|
| chr5 *FGF10* TAD | 612 / 34 (block of 607 variants, 43,600,008–44,824,921) | 265 / 1,493 | 544 / 1,656 (run did not finish) |
| chr10 *FGFR2* | 125 / 9 (block of 124 variants, 121,482,293–121,598,063) | 69 / 202 | 102 / 208 |
| chr17 *TBX4* | 273 / 10 (block of 257 variants, 61,201,667–61,484,743) | 126 / 430 | 229 / 322 |

Reads covering the phased variants in AD143 rose from 207 to 1,325 (chr10), 455 to 3,254 (chr17) and
634 to 5,926 (chr5) between the two 6-sample runs, consistent with the resequenced library. The WhatsHap windows for chr10 and chr17 (gene / super-enhancer to gene) are narrower than the
windows used for the site-level inheritance analysis (`config/loci.bed`), which come from the joint VCF.

**Crossovers.** Discriminating non-damage sites are ordered by position; a switch between long runs of
concordant and discordant sites indicates paternal recombination. Isolated single-site switches are
treated as genotyping error.

**Haplotype labels** are locus-specific: "A" is the paternal haplotype carried by AD144.

**Annotation.** VEP 115.1 (offline cache, one consequence per variant), gnomAD MAX_AF (rare: < 0.01),
CADD v1.7 PHRED (≥ 10 highlighted), UCSC phyloP100way (≥ 1.0 = conserved). For indels, phyloP is the
maximum over the reference bases spanned; indels absent from the CADD gnomAD indel file are reported as NA.


## Analyses performed outside this code

- **AlphaGenome** variant prioritization (AlphaGenome Variant Impact score) — [(https://deepmind.google/science/alphagenome/)].
- **FABIAN-variant** (2026) transcription-factor binding predictions — [web interface](https://fabianapp.org/variant26/).
- **Hi-C** (IMR-90, GEO GSE63525) and ENCODE ChIP-seq tracks — visualized in [[browser](https://3dgenome.fsm.northwestern.edu/)]; no processing code.

## Reproducibility notes

These scripts are cleaned-up versions of the commands run interactively during the analysis.
Before release, rerun them on the original data and confirm that they reproduce the values above.

- `06_fgf10_hapA_candidates.sh` reads the joint VCF directly; the original selection was made from a
  genotype matrix derived from the same VCF. Confirm that it returns the same 30 sites.
- `08_vep.sh` uses `--pick` with an explicit `--pick_order`. The original *FGFR2* run used `--pick`
  with VEP's default order, and the original *FGF10* run used `--flag_pick_allele_gene` with the order
  shown here. Allele frequencies are unaffected; the selected consequence can differ for a few variants.
- `01_joint_genotyping.sh` reproduces the per-chromosome joint-calling run; per-sample gVCFs came from
  bcbio-nextgen (HaplotypeCaller, GVCF mode).

## Citation

[Citation / CITATION.cff]

## License

[Choose a license, e.g. MIT, and add a LICENSE file.]
