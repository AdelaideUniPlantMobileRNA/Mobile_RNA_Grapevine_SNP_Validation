#!/bin/bash
#SBATCH --job-name=master_db_v3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --partition=batch
# Step 2.6 (v3) wrapper: 09 v3 -> results/dipcall/master_diagnostic_db_v3.tsv
# (donor allele confirmed in ANY donor tissue). The v1 database is not touched.
set -euo pipefail
module load Anaconda3/2025.06-1
source activate /scratchdata1/users/a1709557/conda_envs/grapevine_snp

PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
GFF=/scratchdata1/users/a1709557/Own_RNA/Reference_genomes/Grape_reference_annotations/PN40024_5.1_on_T2T_ref_with_names.gff3

Rscript "$PROJ/scripts/09_build_master_db_v3.R" \
    "$PROJ/results/dipcall" \
    "$PROJ/results/rna/cnt" \
    "$PROJ/metadata/step2_samples.tsv" \
    "$GFF" \
    10 3 0.10 2 10 0.02
