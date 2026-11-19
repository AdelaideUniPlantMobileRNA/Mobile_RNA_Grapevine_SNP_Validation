#!/bin/bash
#SBATCH --job-name=rna_discover
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --array=1-28
# =====================================================================
# Step 2.3 : RNA variant DISCOVERY per sample (homograft).
# bcftools mpileup over the exons, keeping only positions where this sample
# has a non-reference allele (ALT != "<*>"). These seed the master position
# pool so that assembly-gap and clone-specific variants (invisible to the
# genome assemblies) are recovered. Full filtering (L1-L4) happens later in 09.
# Output <sample>.disc.tsv:  CHROM  POS  REF  ALT  AD  ADF  ADR
# Submit after 05 (exons.bed) and 04 (BAMs).
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
EXONS=$RNA/exons.bed
DISC=$RNA/disc
mkdir -p "$DISC"

LINE=$(( SLURM_ARRAY_TASK_ID + 1 ))
read -r SAMPLE GROUP TISSUE REP R1 R2 < <(sed -n "${LINE}p" "$SHEET")
echo ">> [$SLURM_ARRAY_TASK_ID] discover variants in $SAMPLE ($GROUP/$TISSUE)"

# Discovery only needs to SEE that an alt allele is present (precise counting is
# done later in 08), so cap depth low (-d 1000) -> far faster at high-coverage genes.
bcftools mpileup -f "$REF" -R "$EXONS" \
     -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR -q 20 -Q 20 -B --max-depth 1000 \
     -Ou "$ALIGN/${SAMPLE}.bam" \
 | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%ADF\t%ADR]\n' \
 | awk -F'\t' '$4!="<*>"' \
 > "$DISC/${SAMPLE}.disc.tsv"

echo ">> wrote $DISC/${SAMPLE}.disc.tsv  ($(wc -l < "$DISC/${SAMPLE}.disc.tsv") alt-bearing positions)"
