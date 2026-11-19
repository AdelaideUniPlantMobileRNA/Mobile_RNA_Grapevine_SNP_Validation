#!/bin/bash
#SBATCH --job-name=star_hetero
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --time=08:00:00
#SBATCH --partition=batch
#SBATCH --array=1-32            # = number of heterograft samples (SRscion12 + RSscion12 + SRR4 + RSR4)
# =====================================================================
# Step 4.1 : STAR-align each heterograft RNA sample to the SAME T2T index
# used for the homografts (common backbone -> same diagnostic coordinates).
# Same STAR settings as 04 (unique mappers only, two-pass, read groups).
# Submit after make_hetero_sample_sheet.sh.
# =====================================================================
set -euo pipefail
module load Anaconda3/2025.06-1
source activate /scratchdata1/users/a1709557/conda_envs/grapevine_snp

PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
INDEX=$PROJ/results/star_index
SHEET=$PROJ/metadata/hetero_samples.tsv
ALIGN=$PROJ/results/star_hetero
THREADS=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$ALIGN"

# columns: sample graft local_geno mobile_geno direction tissue rep R1 R2
LINE=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r SAMPLE GRAFT LG MG DIR TISSUE REP R1 R2 < <(sed -n "${LINE}p" "$SHEET")
echo ">> [$SLURM_ARRAY_TASK_ID] $SAMPLE ($GRAFT local=$LG detect=$MG/$DIR $TISSUE rep$REP)"

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
