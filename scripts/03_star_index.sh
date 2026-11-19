#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --partition=batch
# =====================================================================
# Step 2.0 : build the STAR genome index on the T2T backbone (one-time).
# Uses the PN40024-on-T2T annotation so splice junctions are known.
# =====================================================================
set -euo pipefail
module load Anaconda3/2025.06-1
source activate /scratchdata1/users/a1709557/conda_envs/grapevine_snp

# ---------- paths ----------
GENODIR=/scratchdata1/users/a1709557/Own_RNA/Reference_genomes/Grape_reference_assemblies
REF=$GENODIR/T2T_ref.fasta
GFF=/scratchdata1/users/a1709557/Own_RNA/Reference_genomes/Grape_reference_annotations/PN40024_5.1_on_T2T_ref_with_names.gff3
PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
INDEX=$PROJ/results/star_index

# ---------- parameters ----------
THREADS=${SLURM_CPUS_PER_TASK:-16}
# sjdbOverhang = (max read length - 1). 2x150 bp reads -> 149. Set to your read length - 1.
SJDB_OVERHANG=149
# genomeSAindexNbases: min(14, log2(genome_len)/2 - 1). ~500 Mb grape genome -> 13.
SA_NBASES=13

mkdir -p "$INDEX"
STAR --runMode genomeGenerate \
     --runThreadN "$THREADS" \
     --genomeDir "$INDEX" \
     --genomeFastaFiles "$REF" \
     --sjdbGTFfile "$GFF" \
     --sjdbGTFtagExonParentTranscript Parent \
     --sjdbGTFfeatureExon exon \
     --sjdbOverhang "$SJDB_OVERHANG" \
     --genomeSAindexNbases "$SA_NBASES"

echo ">> STAR index built in $INDEX"
