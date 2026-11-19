#!/bin/bash
#SBATCH --job-name=detect_v3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --partition=batch
#SBATCH --array=1-32            # one task per heterograft sample
# =====================================================================
# Step 4a (v3) wrapper: per-SNP AND gene-pooled baymobil Bayes factors,
# one heterograft sample per task. Self-graft cap 2000, SNP thinning 300 bp.
# Needs results/dipcall/master_diagnostic_db_v3.tsv (run 10_run_master_db_v3.sh first).
# Results: results/mobile_v3 (v1, v2 kept). Finished samples are skipped.
# =====================================================================
set -euo pipefail
ENVPY=/scratchdata1/users/a1709557/conda_envs/grapevine_snp/bin/python
PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
SHEET=$PROJ/metadata/hetero_samples.tsv

SAMPLE=$(sed -n "$(( SLURM_ARRAY_TASK_ID + 1 ))p" "$SHEET" | cut -f1)

"$ENVPY" "$PROJ/scripts/13_detect_mobile_v3.py" \
    "$PROJ/results/dipcall/master_diagnostic_db_v3.tsv" \
    "$PROJ/results/dipcall/diagnostic_sites.tsv" \
    "$PROJ/results/rna/cnt" \
    "$PROJ/results/rna/hetero_cnt" \
    "$PROJ/metadata/step2_samples.tsv" \
    "$SHEET" \
    "$PROJ/results/mobile_v3" \
    "$SAMPLE" \
    10 2000 300
