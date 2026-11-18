# ---------------------------------------------------------------------------
# config.sh — edit these paths for your system, then every script sources it.
# ---------------------------------------------------------------------------

# Reference genome (GRCh38, indexed .fa or .fa.gz with .fai/.gzi)
REF=/path/to/Homo_sapiens.GRCh38.dna.toplevel.fa.gz

# GATK launcher and threads
GATK=gatk
THREADS=4

# Joint VCF produced by 01_joint_genotyping.sh (or an existing one)
JOINT_VCF=/path/to/family.joint_v2.vcf.gz

# BAM files (used for coverage and WhatsHap read-based phasing)
#   BAM_TEMPLATE: {S} is replaced by the sample ID
BAM_TEMPLATE=/path/to/final/{S}_sort_bam/{S}_sort_bam-ready.bam
AD143_BAM_ORIGINAL=/path/to/original/AD143_sort_bam-ready.bam
AD143_BAM_RESEQ=/path/to/resequenced/AD143_3-ready.bam

# Sample names in the VCF/BAM read groups are "<ID><SUFFIX>", e.g. AD141_sort_bam
SAMPLE_SUFFIX=_sort_bam

# Pedigree roles
FATHER=AD141
MOTHER=AD142
PROBAND=AD143
SIBS="AD144 AD145 AD146"
AFFECTED_AUNTS="AD137,AD138,AD139,AD140"           # affected relatives outside the sibship
AFFECTED_ALL="AD137,AD138,AD139,AD140,AD142,AD146"  # all affected relatives except the proband

# Ensembl VEP (Docker)
VEP_IMAGE=ensemblorg/ensembl-vep:release_115.1
VEP_CACHE=/path/to/vep_data                 # contains homo_sapiens/115_GRCh38
VEP_FASTA_IN_CACHE=homo_sapiens/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
VEP_PLUGINS=/path/to/vep_data/Plugins      # optional
DBNSFP_DIR=/path/to/plugins_data           # optional; contains dbNSFP4.1a_variant.chr<N>.bgz

# CADD v1.7 precomputed scores (GRCh38; chromosomes without "chr")
CADD_SNV=/path/to/CADD/whole_genome_SNVs.tsv.gz
CADD_INDEL=/path/to/CADD/gnomad.genomes.r4.0.indel.tsv.gz

# UCSC phyloP bigWigs (hg38)
PHYLOP100=/path/to/hg38.phyloP100way.bw
PHYLOP30=/path/to/hg38.phyloP30way.bw

# Output folder
OUT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/results
mkdir -p "$OUT"
export CADD_SNV CADD_INDEL PHYLOP100 PHYLOP30
