#!/bin/bash
#SBATCH --job-name=master_positions
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --partition=batch
# =====================================================================
# Step 2.4 : build the MASTER position list to count at.
#   union of:
#     (a) genome exonic diagnostic positions (Step 1)
#     (b) RNA-discovered positions seen in >= 2 samples (coarse seed; the
#         strict L1-L4 filtering is applied later in 09)
# This is a loose seed on purpose — 09 does the real QC — but requiring 2
# samples here already drops most single-library noise before the count pass.
# Run after 06 (all disc.tsv).
# =====================================================================
set -euo pipefail

PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
DIAG=$PROJ/results/dipcall/diagnostic_sites.tsv
RNA=$PROJ/results/rna
DISC=$RNA/disc
OUT=$RNA/master.regions.tsv

# Keep sort's temp files on scratch, NOT the tiny node-local /tmp (which is what
# ran out of space). -S gives a big in-memory buffer so fewer temp files.
TMP=$RNA/sorttmp; mkdir -p "$TMP"; export TMPDIR="$TMP"
SORTOPT="-T $TMP -S 6G --parallel=${SLURM_CPUS_PER_TASK:-4}"

# (a) genome exonic diagnostic positions (col 14 = is_exonic)
awk -F'\t' 'NR>1 && $14=="TRUE"{print $1"\t"$2}' "$DIAG" > "$RNA/_genome.pos"

# (b) RNA positions with a QUALITY alt call in >= 2 samples.
#     A raw disc.tsv row only means ">=1 non-ref read" (mostly single-read
#     sequencing error). Here we require, per sample, a real alt: site depth>=10,
#     alt reads>=3, VAF>=0.10, and both strands (when alt>=4) -- the same L1 gate
#     09 uses. disc.tsv already carries AD/ADF/ADR, so no re-run of 06 is needed.
#     This collapses the noisy ~17M seed to the positions worth counting.
awk -F'\t' '{
  nalt=split($4,alt,","); split($5,ad,","); split($6,adf,","); split($7,adr,",");
  dp=0; for(i=1;i<=length(ad);i++) dp+=ad[i];
  if(dp<10) next;
  for(k=1;k<=nalt;k++){
    if(alt[k]=="<*>") continue;
    c=ad[k+1]; f=adf[k+1]; r=adr[k+1];
    if(c>=3 && c/dp>=0.10 && (c<4 || (f>=1 && r>=1))){ print $1"\t"$2; break }
  }
}' "$DISC"/*.disc.tsv | sort $SORTOPT | uniq -c \
  | awk '$1>=2{print $2"\t"$3}' > "$RNA/_rna.pos"

cat "$RNA/_genome.pos" "$RNA/_rna.pos" | sort $SORTOPT -k1,1 -k2,2n -u > "$OUT"
rm -f "$RNA/_genome.pos" "$RNA/_rna.pos"
rm -rf "$TMP"

echo ">> master positions: $(wc -l < "$OUT")"
