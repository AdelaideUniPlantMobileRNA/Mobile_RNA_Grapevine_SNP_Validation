#!/bin/bash
#SBATCH --job-name=dipcall_SR
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --output=logs/dipcall_%j.log
#SBATCH --error=logs/dipcall_%j.err
# =====================================================================
# Step 1: build the genome-based diagnostic catalogue.
# dipcall(Li 2018) -> diploid VCF + confident BED of each genotype vs T2T;
# then co-genotype (merge) on T2T coordinates for 02_build_diagnostic_sites.R
# =====================================================================
set -euo pipefail

# ---------- environment ----------
module load Anaconda3/2025.06-1
source activate /scratchdata1/users/a1709557/conda_envs/grapevine_snp

# ---------- paths ----------
GENODIR=/scratchdata1/users/a1709557/Own_RNA/Reference_genomes/Grape_reference_assemblies
T2T_REF=$GENODIR/T2T_ref.fasta
SHZ_PRIMARY=$GENODIR/SHZ_1654_primary.fasta
SHZ_HAP=$GENODIR/SHZ_1654_haplotigs.fasta
RUG_PRIMARY=$GENODIR/Primary_contigs_140Ruggeri.fasta
RUG_ALT=$GENODIR/alt_contigs_140Ruggeri.fasta

PROJ=/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall
OUTDIR=$PROJ/results/dipcall
THREADS=${SLURM_CPUS_PER_TASK:-16}
PER=$(( THREADS / 2 ))          # 2 haplotype alignments in parallel (make -j2)

mkdir -p "$OUTDIR" "$PROJ/logs"
cd "$PROJ"

samtools faidx "$T2T_REF"

# ---------- Shiraz (V. vinifera): default asm5 ----------
run-dipcall -t "$PER" "$OUTDIR/shiraz"  "$T2T_REF" "$SHZ_PRIMARY" "$SHZ_HAP" > "$OUTDIR/shiraz.mak"
make -j2 -f "$OUTDIR/shiraz.mak"
#   -> $OUTDIR/shiraz.dip.vcf.gz + $OUTDIR/shiraz.dip.bed

# ---------- Ruggeri (interspecific): relax asm5 -> asm20 ----------
run-dipcall -t "$PER" "$OUTDIR/ruggeri" "$T2T_REF" "$RUG_PRIMARY" "$RUG_ALT" > "$OUTDIR/ruggeri.mak"
sed -i 's/-x asm5/-x asm20/g; s/-xasm5/-xasm20/g' "$OUTDIR/ruggeri.mak"
echo ">> Ruggeri preset check (want asm20):"; grep -o 'asm[0-9]*' "$OUTDIR/ruggeri.mak" | sort -u
make -j2 -f "$OUTDIR/ruggeri.mak"
#   -> $OUTDIR/ruggeri.dip.vcf.gz + $OUTDIR/ruggeri.dip.bed

# ---------- PRIMARY (hap1) callable BED for each genotype ----------
# dipcall's *.dip.bed needs BOTH haplotypes aligned, so it discards the
# homozygous majority of a Purge-Haplotigs primary assembly. 02 instead
# uses the PRIMARY (hap1) callable region: the "R" (region) records that
# paftools emits into *.hap1.var.gz are exactly those confident regions.
for g in shiraz ruggeri; do
  if [[ -s "$OUTDIR/${g}.hap1.bed" ]]; then
    echo ">> ${g}.hap1.bed already present (dipcall)."
  else
    gzip -dc "$OUTDIR/${g}.hap1.var.gz" \
      | awk 'BEGIN{OFS="\t"} $1=="R"{print $2,$3,$4}' \
      | sort -k1,1 -k2,2n > "$OUTDIR/${g}.hap1.bed"
  fi
  echo "   ${g}.hap1.bed bp: $(awk '{s+=$3-$2} END{print s}' "$OUTDIR/${g}.hap1.bed")"
done

# ---------- reheader samples, normalise, index ----------
for g in shiraz ruggeri; do
  echo "$g" > "$OUTDIR/${g}.sample"
  bcftools reheader -s "$OUTDIR/${g}.sample" "$OUTDIR/${g}.dip.vcf.gz" \
    | bcftools norm -f "$T2T_REF" -m -any -Oz -o "$OUTDIR/${g}.norm.vcf.gz"
  bcftools index -t "$OUTDIR/${g}.norm.vcf.gz"
done

# ---------- merge on T2T coordinates (missing stays ./.) ----------
bcftools merge -Oz -o "$OUTDIR/merged.vcf.gz" \
  "$OUTDIR/shiraz.norm.vcf.gz" "$OUTDIR/ruggeri.norm.vcf.gz"
bcftools index -t "$OUTDIR/merged.vcf.gz"

# ---------- flat genotype table for R ----------
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$OUTDIR/merged.vcf.gz" \
  > "$OUTDIR/merged_genotypes.tsv"
# columns: CHROM POS REF ALT shiraz_GT ruggeri_GT

echo ">> DONE. Key outputs in $OUTDIR :"
echo "   shiraz.dip.bed  ruggeri.dip.bed  merged_genotypes.tsv"
