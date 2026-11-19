# Steps 2–5 — RNA markers, Bayesian detection, verification

All scripts are SLURM jobs for Phoenix (paths inside the scripts point to
`/scratchdata1/users/a1709557/Own_RNA/Validation_withSNP_Calling_dipcall`).
Every threshold, with its source (published vs our choice), is listed in
`docs/assumptions.md`. The history of changes and errors is in `docs/analysis_log.md`.

## Run order (final version)

```bash
# ---- Step 1 (genome markers) -------------------------------------------------
sbatch scripts/01_dipcall.sh
Rscript scripts/02_build_diagnostic_sites.R results/dipcall <PN40024_on_T2T.gff3>
#   -> results/dipcall/diagnostic_sites.tsv   (test: bash scripts/tests/run_test.sh)

# ---- Step 2 (homograft RNA) --------------------------------------------------
sbatch scripts/03_star_index.sh
bash   scripts/make_step2_sample_sheet.sh          # metadata/step2_samples.tsv (28 samples)
sbatch scripts/04_star_align.sh                    # array 1-28, STAR two-pass, unique reads
sbatch scripts/05_make_exons_bed.sh                # exons.bed from the GFF3
sbatch scripts/06_rna_discover.sh                  # array 1-28, alt-bearing exonic positions
sbatch scripts/07_build_positions.sh               # master.regions = genome exonic ∪ RNA(>=2 samples)
sbatch scripts/08_rna_count.sh                     # array 1-28, AD/ADF/ADR at master positions
sbatch scripts/10_run_master_db_v3.sh              # runs 09_build_master_db_v3.R
#   -> results/dipcall/master_diagnostic_db_v3.tsv (test: bash scripts/tests/run_test_09_v3.sh)

# ---- Step 4 (heterografts, Bayesian) -----------------------------------------
bash   scripts/make_hetero_sample_sheet.sh         # metadata/hetero_samples.tsv (32 samples)
sbatch scripts/11_star_align_hetero.sh             # array 1-32
sbatch scripts/12_hetero_count.sh                  # array 1-32
sbatch scripts/14_run_detect_v3.sh                 # array 1-32, runs 13_detect_mobile_v3.py
sbatch scripts/15_aggregate_mobile_v3.sh           # 13b_aggregate_mobile_v3.py + 16_inspect_mobile.py

# ---- Step 5 (verification) ---------------------------------------------------
python scripts/17_replicate_concordance.py results/mobile_v3
python scripts/18_contamination_check.py results/mobile_v3 results/rna/cnt results/rna/hetero_cnt \
       metadata/step2_samples.tsv metadata/hetero_samples.tsv SRR4      # run with sbatch, 64G
python scripts/19_spotcheck.py SRL1 20 1          # samtools re-count + baymobil run_bayes
python scripts/check_gene_labels.py <gff3> results/dipcall/diagnostic_sites.tsv \
       results/dipcall/master_diagnostic_db_v3.tsv
```
Sample SRR4 is excluded from all heterograft summaries (label error; see analysis_log).

## What each Step-2/4/5 script does

| Script | Input → output | Key point |
|---|---|---|
| 06 | homograft BAM → `results/rna/disc/*.disc.tsv` | bcftools mpileup over exons, alt-bearing positions only (discovery) |
| 07 | diagnostic_sites + disc → `master.regions.tsv` | union of genome exonic markers and RNA positions passing the L1 gate in ≥2 samples |
| 08 / 12 | BAM → `results/rna/cnt/` / `hetero_cnt/` | allele **counting** at known positions (`bcftools mpileup -T`, `-q20 -Q20 -B -d10000`), no genotype calling |
| 09 v3 | homograft counts → `master_diagnostic_db_v3.tsv` | hard filters L1–L4 (below); donor allele from **any** donor tissue |
| 13 v3 | markers + counts → `per_sample/<s>.snp.tsv`, `<s>.gene.tsv` | baymobil `fasterpostN2` per SNP; gene-pooled BF on SNPs ≥300 bp apart |
| 13b v3 | per-sample files → `mobile_transcripts*.tsv` | per-SNP sum (published) and gene-pooled (ours), side by side |
| 16 | → `check_global_fraction.tsv`, `check_mobile_genes.tsv` | global donor fraction het vs hom; per-gene detail of calls |
| 17 | → `check_concordance_*.tsv` | observed vs chance-expected replicate agreement |
| 18 | → `check_contamination*.tsv` | n ~ Poisson(a·N + c·D); heterograft vs homograft background |
| 19 | → `spotcheck_<sample>.tsv` | independent samtools counts + baymobil `run_bayes` |
| check_gene_labels | any table with CHROM, POS, gene_id | every gene_id must contain its position (GFF3) |

## Homograft filters (09 v3), per genotype SS / RR

| Level | Rule |
|---|---|
| L1 per sample | site depth ≥ 10, allele reads ≥ 3, allele fraction ≥ 0.10; if allele reads ≥ 4, seen on both read orientations (ADF and ADR ≥ 1) |
| L2 reproducibility | passes L1 in ≥ 2 replicates of the same genotype × tissue, and both orientations in the pooled reads |
| L3 RNA editing | drop A→G on (+)-strand genes and T→C on (−)-strand genes when the allele has no genome support |
| L4 donor-private | donor allele passes L1–L3 in **any** donor tissue; recipient same tissue: pooled depth ≥ 10, ≥ 2 replicates covered, allele fraction ≤ 2% |

## Detection rules (13 v3 + 13b v3)

* Marker tested in a heterograft sample if it has ≥ 1 read there (no homograft
  expression requirement).
* Marker removed if the tissue-matched recipient homograft shows the donor allele
  (≥ 1 read and > 2% of reads), at any depth.
* baymobil inputs: N, n = heterograft total / donor-allele reads; Nh1, nh1 = recipient
  homograft (same tissue, replicates pooled) total / donor-allele reads; Nh2, nh2 =
  donor homograft (all tissues pooled) total / non-donor reads. Homograft depth capped
  at 2000 (counts scaled). nmax = 10. No homograft reads → flat Beta(1,1) prior.
* Gene log10BF in a replicate = sum of its SNP log10BF; replicate passes at ≥ 1;
  mobile if ≥ 2 replicates pass.

## Main output files (`results/mobile_v3/`, not in git)

| File | Content |
|---|---|
| `mobile_transcripts.tsv` | every gene × graft × tissue × direction tested; `mobile` = True for calls |
| `mobile_transcripts_per_rep.tsv` | per gene × sample: n_snp, sum_log10BF, pass |
| `mobile_transcripts_both.tsv` | per-SNP sum and gene-pooled calls side by side |
| `per_sample/<sample>.snp.tsv` | per SNP: position, donor allele, N, n, Nh1, nh1, Nh2, nh2, log10BF |
| `candidates_contamination_robust.tsv` | high-confidence subset (obs/exp ≥ 5 for both sources), wet-lab priority only |
