#!/bin/bash
# =====================================================================
# Environment setup on an HPC with the "modules + Anaconda3" layout
# (tested on the University of Adelaide Phoenix HPC).
#
# Builds one conda environment holding every tool used by the pipeline:
#   dipcall, minimap2, k8            (Step 1.1 genome diagnostic calling)
#   samtools, bcftools, bedtools     (variant handling)
#   r-base, r-data.table             (Step 1.2 diagnostic-site builder)
#   star, gatk4                      (Step 2 RNA confirmation - later)
#
# Run ONCE, interactively, on a login node (it only downloads/installs):
#   bash 00_setup_env.sh
# =====================================================================
set -euo pipefail

# ---------------------------------------------------------------------
# EDIT THESE THREE PATHS for your account, then run.
# ---------------------------------------------------------------------
ENV_PREFIX=/scratchdata1/users/a1709557/conda_envs/grapevine_snp   # where the env lives
PKG_CACHE=/scratchdata1/users/a1709557/conda_pkgs                  # package cache (keep off /home)
ANACONDA_MODULE=Anaconda3/2025.06-1                                # `module avail Anaconda3` to check
# ---------------------------------------------------------------------

# Why a custom prefix + cache:
#   * /home has a small quota; conda envs and the package cache are large,
#     so both are placed on /scratchdata1.
#   * `module load` gives us the base conda; we never install into it.
module load "$ANACONDA_MODULE"
export CONDA_PKGS_DIRS="$PKG_CACHE"
mkdir -p "$PKG_CACHE"

# Why --override-channels -c conda-forge -c bioconda:
#   the shared Anaconda "defaults" channel triggers a non-interactive
#   Terms-of-Service error on managed installs; overriding the channels
#   avoids `defaults` entirely and is also the recommended bioconda setup.
conda create -y -p "$ENV_PREFIX" --override-channels \
  -c conda-forge -c bioconda \
  dipcall minimap2 k8 \
  samtools bcftools bedtools \
  r-base r-data.table \
  star gatk4

# Activate by PATH (the env is not registered under a name).
source activate "$ENV_PREFIX"

echo "===================================================================="
echo "Verifying tools:"
for t in run-dipcall minimap2 k8 samtools bcftools bedtools STAR gatk Rscript; do
  printf '  %-12s ' "$t"; command -v "$t" || { echo "MISSING"; exit 1; }
done
Rscript -e 'library(data.table); cat("  data.table loads OK:", as.character(packageVersion("data.table")), "\n")'
echo "===================================================================="
echo "Done. In every job / session, activate with:"
echo "    module load $ANACONDA_MODULE"
echo "    source activate $ENV_PREFIX"
