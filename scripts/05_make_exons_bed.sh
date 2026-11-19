#!/bin/bash
#SBATCH --job-name=exons_bed
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --partition=batch
# =====================================================================
# Step 2.2 : extract the exon intervals from the T2T-backbone annotation.
# RNA variant discovery (06) is restricted to these (mobile mRNA is exonic,
# and it keeps the mpileup tractable). Run once.
# =====================================================================
set -euo pipefail
module load Anaconda3/2025.06-1
source activate /scratchdata1/users/a1709557/conda_envs/grapevine_snp

GFF=/scratchdata1/users/a1709557/Own_RNA/Reference_genomes/Grape_reference_annotations/PN40024_5.1_on_T2T_ref_with_names.gff3
PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
OUT=$PROJ/results/rna
mkdir -p "$OUT"

awk -F'\t' '$3=="exon"{print $1"\t"($4-1)"\t"$5}' "$GFF" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i - > "$OUT/exons.bed"

echo ">> exons: $(wc -l < "$OUT/exons.bed") merged intervals, $(awk '{s+=$3-$2} END{print s}' "$OUT/exons.bed") bp"
