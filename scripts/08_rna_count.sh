#!/bin/bash
#SBATCH --job-name=rna_count
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --array=1-28
# =====================================================================
# Step 2.5 : count per-allele reads (with strand) at the MASTER positions,
# for every homograft sample. Unlike discovery (06) this is NOT alt-filtered,
# so reference-only positions are included too — needed so 09 knows the
# recipient's depth / background where it carries no alternate allele.
# Output <sample>.cnt.tsv:  CHROM  POS  REF  ALT  AD  ADF  ADR
# Submit after 07 (master.regions.tsv).
# =====================================================================
set -euo pipefail
module load Anaconda3/2025.06-1
source activate /scratchdata1/users/a1709557/conda_envs/grapevine_snp

GENODIR=/scratchdata1/users/a1709557/Own_RNA/Reference_genomes/Grape_reference_assemblies
REF=$GENODIR/T2T_ref.fasta
PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
SHEET=$PROJ/metadata/step2_samples.tsv
ALIGN=$PROJ/results/star
RNA=$PROJ/results/rna
REGIONS=$RNA/master.regions.tsv
CNT=$RNA/cnt
mkdir -p "$CNT"

LINE=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r SAMPLE GROUP TISSUE REP R1 R2 < <(sed -n "${LINE}p" "$SHEET")
echo ">> [$SLURM_ARRAY_TASK_ID] count at master positions for $SAMPLE"

# -T (targets, STREAM through the BAM) not -R (regions, index-seek per position):
# with ~1.45M scattered single positions, -R does ~1.45M index seeks and is
# pathologically slow; -T scans the BAM once -> orders of magnitude faster,
# identical counts. -d 10000 is ample for an allele fraction.
bcftools mpileup -f "$REF" -T "$REGIONS" \
     -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR -q 20 -Q 20 -B --max-depth 10000 \
     -Ou "$ALIGN/${SAMPLE}.bam" \
 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%ADF\t%ADR]\n' \
 > "$CNT/${SAMPLE}.cnt.tsv"

echo ">> wrote $CNT/${SAMPLE}.cnt.tsv  ($(wc -l < "$CNT/${SAMPLE}.cnt.tsv") positions)"
