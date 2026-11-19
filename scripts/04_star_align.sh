#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --partition=batch
#SBATCH --array=1-28            # = number of SS + RR samples (16 + 12)
# =====================================================================
# Step 2.1 : STAR align each homograft RNA sample to the T2T index.
#   * unique mappers only (ASE requires it): --outFilterMultimapNmax 1
#   * two-pass for better splice-junction detection
#   * read-group tag per sample (needed by GATK downstream)
# Submit AFTER make_step2_sample_sheet.sh has written metadata/step2_samples.tsv.
#   sbatch 04_star_align.sh
# =====================================================================
set -euo pipefail
module load Anaconda3/2025.06-1
source activate /scratchdata1/users/a1709557/conda_envs/grapevine_snp

PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
INDEX=$PROJ/results/star_index
SHEET=$PROJ/metadata/step2_samples.tsv       # or repo metadata/ copied in
ALIGN=$PROJ/results/star
THREADS=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$ALIGN"

# ---- pick this array task's sample (skip header) ----
LINE=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r SAMPLE GROUP TISSUE REP R1 R2 < <(sed -n "${LINE}p" "$SHEET")
echo ">> [$SLURM_ARRAY_TASK_ID] $SAMPLE ($GROUP/$TISSUE rep$REP)"
echo "   R1=$R1"; echo "   R2=$R2"

STAR --runThreadN "$THREADS" \
     --genomeDir "$INDEX" \
     --readFilesIn "$R1" "$R2" \
     --readFilesCommand zcat \
     --twopassMode Basic \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS nM MD \
     --outSAMattrRGline ID:"$SAMPLE" SM:"$SAMPLE" PL:ILLUMINA LB:"$SAMPLE" \
     --outFileNamePrefix "$ALIGN/${SAMPLE}."

mv "$ALIGN/${SAMPLE}.Aligned.sortedByCoord.out.bam" "$ALIGN/${SAMPLE}.bam"
samtools index "$ALIGN/${SAMPLE}.bam"
echo ">> done: $ALIGN/${SAMPLE}.bam"
