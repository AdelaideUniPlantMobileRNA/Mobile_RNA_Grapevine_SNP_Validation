#!/bin/bash
#SBATCH --job-name=agg_v3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --partition=batch
# Step 4b+4c (v3): two summaries side by side (per-SNP sum, gene-pooled), then the
# quality checks of 16 on the per-SNP-sum calls.
set -euo pipefail
ENVPY=/scratchdata1/users/a1709557/conda_envs/grapevine_snp/bin/python
PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
OUT=$PROJ/results/mobile_v3
NS=$(ls "$OUT"/per_sample/*.snp.tsv 2>/dev/null | wc -l)
NG=$(ls "$OUT"/per_sample/*.gene.tsv 2>/dev/null | wc -l)
echo ">> per-sample files: snp $NS / 32, gene $NG / 32"
[[ "$NS" -eq 32 && "$NG" -eq 32 ]] || { echo "Not all samples finished - resubmit 14_run_detect_v3.sh first."; exit 1; }
EXCLUDE="SRR4"   # suspected label error (labelled Ruggeri root, mostly Shiraz RNA)
"$ENVPY" "$PROJ/scripts/13b_aggregate_mobile_v3.py" "$OUT" 2 "$EXCLUDE"
echo; echo "################ quality checks (per-SNP-sum calls) ################"
"$ENVPY" "$PROJ/scripts/16_inspect_mobile.py" "$OUT" "$PROJ/results/mobile_v2/mobile_transcripts.tsv" "$EXCLUDE" 2000
