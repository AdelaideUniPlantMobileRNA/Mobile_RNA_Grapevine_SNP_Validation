# Analysis log

> A dated, append-only lab notebook for the computational work.
> One entry per working session. Record: what you ran, on what input, which script,
> the key result, and any decision. This is what saves you when a reviewer asks
> "how exactly did you get this number?".

---

### 2026-10-20  — set up project + environment
- Created project structure, git repo, conda env (environment.yml).
- Commit: `a1b2c3d`
- Decision: using mamba; tool versions pinned.

### YYYY-MM-DD  — <short title>
- **Goal:**
- **Input:** path/to/input (describe)
- **Command / script:** `scripts/0X_xxx.sh`  (SLURM job id: 123456)
- **Parameters that matter:**
- **Output:** path/to/output
- **Result / observation:**
- **Decision / next step:**
- **Commit:** `<git hash>`

---

<!-- template to copy for each session:

### YYYY-MM-DD  — title
- Goal:
- Input:
- Script/command:
- Output:
- Result:
- Decision:
- Commit:

-->

---
## Log entries (Oct–Nov 2026)

### 2026-10-28 — Step 1 genome markers (dipcall + 02)
- dipcall Shiraz (asm5) / Ruggeri (asm20) vs T2T; merged genotypes; 02 builds allele sets and private alleles.
- **Later found to be affected by a row-order bug (see 2026-11-16); results of this date are superseded.**

### 2026-11-05 to 11-09 — Step 2 homograft RNA (03–09)
- STAR two-pass, unique reads; 28 homograft samples (SS/RR × B, L, PY, R).
- GATK ASEReadCounter dropped ("More then one variant context" error) → bcftools mpileup AD counts.
- 07: an L1 gate was added (16.9M → 1.45M positions); `sort -T` moved to scratch (node /tmp full).
- 08: switched `-R` to `-T` (streaming; `-R` was pathologically slow).

### 2026-11-10 — Step 4 v1 (baymobil)
- 13 called baymobil's `fasterpostN2` directly with a process pool (run_bayes' pandarallel used all 72 cores, >200 GB RAM).
- v1 tested only markers expressed in BOTH homografts of the same tissue (~2000–3600 genes per combination).
- **Superseded** (design error + 2026-11-16 bug).

### 2026-11-11 — sample SRR4 excluded
- SRR4 (labelled SR graft, Ruggeri root) carried mostly Shiraz RNA (distal fraction ≈ 0.71) → label error / mix-up; excluded from all heterograft summaries. To be re-confirmed with corrected markers.

### 2026-11-11 to 11-13 — Step 4 v2
- Test every marker with ≥1 heterograft read (no homograft-expression requirement); genome markers taken directly from Step 1.
- Fixed double counting of depth at positions with an extra mpileup indel record.
- **Superseded** (2026-11-16 bug).

### 2026-11-13 — Step 4 v3
- 09 v3: donor allele confirmed in ANY donor tissue.
- Recipient-homograft removal at any depth (v2 required depth ≥10; low-depth sites leaked).
- Homograft cap 300 → 2000; gene-pooled BF added as a side-by-side summary.

### 2026-11-16 — ROW-ORDER BUG found and fixed in 02 and 09
- **Symptom:** a gene (chr17g08320) had 176 markers in v2 and 0 in v3; check against GFF3 showed gene_id matched the position for 0.00% of rows in `diagnostic_sites.tsv`.
- **Cause:** `setkey()` re-sorts the query table before `foverlaps()`; results came back in (chrom, pos) order but were assigned in the original (VCF contig) order, where chr00 is last. Every answer belonged to another row.
- **Affected in 02:** (1) hap1-coverage test used to infer 0/0 for no-calls → the diagnostic site set itself; (2) gene_id; (3) is_exonic (→ 07 positions, 13 marker set). **In 09:** gene_id of RNA-only markers.
- **Fix:** carry a row id, restore order, `stopifnot()` check; tests now include a second chromosome listed first (old code fails).
- **Consequence:** all earlier results (v1 4 genes, v2 40, v3 68) withdrawn; reran 02, 07, 08, 12, 09 v3, 13 v3, 15.
- Corrected 02: 4,211,623 candidate sites (3,610,286 SNPs, 601,337 indels); 549,245 exonic. Gene labels 100% correct (check_gene_labels.py). 331,045 counting positions were new.

### 2026-11-17 — Results (corrected v3)
- 177,421 gene × tissue × direction tested; per-SNP-sum calls (≥2 reps): **396** (pooled: 903; both: 176).
- Upward: SR (140-Ru → Shiraz scion B/L/PY) 122; RS (Shiraz → 140-Ru scion) 190. Downward: SR R 40; RS R 44. ≥3 reps: 24 / 51 / 5 / 14.
- Global donor-allele fraction het/hom: RS 0.66–1.24, SR 1.06–1.40.

### 2026-11-18 — Verification
- Replicate concordance (17): ≥2 reps 396 observed vs 28.8 expected by chance; ≥3: 94 vs 0.1; all 4: 29 vs 0.
- Contamination (18): heterograft share of expression-proportional donor reads not above homograft background → no heterograft-specific contamination. 168 calls with obs/exp ≥5 for both sources = high-confidence subset (prioritisation only).
- Spot check (19): 40 random markers (SRL1, RSR4): samtools counts identical for 39/40 (one >8000× site differs 5% from depth capping); Bayes factors identical to baymobil `run_bayes` 40/40; calls agree 40/40.

### 2026-11-19 — Indels evaluated, not used
- Exonic indels ≤15 bp would add ~14% markers but only 212 (Ruggeri→scion) and 242 (Shiraz→root) genes without SNP markers (~2%). Not implemented.

### Known limitation (documented)
- Per-SNP summation: a covered SNP with n = 0 scores ≈ −3.1 regardless of depth (cap 2000), so genes with many SNPs and low mobile fraction are penalised → results are conservative. Gene-pooled reported as sensitivity analysis.
