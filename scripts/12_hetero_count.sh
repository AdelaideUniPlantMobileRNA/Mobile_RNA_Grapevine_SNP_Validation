#!/bin/bash
#SBATCH --job-name=hetero_count
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --array=1-32
# =====================================================================
# Step 4.2 : count per-allele reads (with strand) at the master diagnostic
# positions for each HETEROGRAFT sample. Same as 08 (-T streaming, -d 10000).
# These give the heterograft (n, k) that Step 4's Bayesian test compares against
# the tissue-matched homograft background.
# Output results/rna/hetero_cnt/<sample>.cnt.tsv : CHROM POS REF ALT AD ADF ADR
# =====================================================================
set -euo pipefail
module load Anaconda3/2025.06-1
source activate /scratchdata1/users/a1709557/conda_envs/grapevine_snp

GENODIR=/scratchdata1/users/a1709557/Own_RNA/Reference_genomes/Grape_reference_assemblies
REF=$GENODIR/T2T_ref.fasta
PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
SHEET=$PROJ/metadata/hetero_samples.tsv
ALIGN=$PROJ/results/star_hetero
RNA=$PROJ/results/rna
REGIONS=$RNA/master.regions.tsv
CNT=$RNA/hetero_cnt
mkdir -p "$CNT"

LINE=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r SAMPLE GRAFT LG MG DIR TISSUE REP R1 R2 < <(sed -n "${LINE}p" "$SHEET")
echo ">> [$SLURM_ARRAY_TASK_ID] count $SAMPLE at master positions"

bcftools mpileup -f "$REF" -T "$REGIONS" \
     -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR -q 20 -Q 20 -B --max-depth 10000 \
     -Ou "$ALIGN/${SAMPLE}.bam" \
 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%ADF\t%ADR]\n' \
 > "$CNT/${SAMPLE}.cnt.tsv"

echo ">> wrote $CNT/${SAMPLE}.cnt.tsv  ($(wc -l < "$CNT/${SAMPLE}.cnt.tsv") positions)"
