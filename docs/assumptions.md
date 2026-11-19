# Assumptions and parameter register

Every choice that affects the result, with its source.
**Published** = taken from a peer-reviewed method; **Ours** = a decision made in this
project (with the reason). Reviewers should treat "Ours" items as open to challenge.

| # | Step / script | Choice | Value | Source | Reason / note |
|---|---|---|---|---|---|
| 1 | 01 | dipcall preset | Shiraz asm5, Ruggeri asm20 | Ours | Ruggeri (non-vinifera hybrid) is more divergent from T2T |
| 2 | 02 | callable region | primary (hap1) coverage, not dip.bed | Ours, based on Onetto et al. 2023 (Purge Haplotigs assemblies) | primary is collapsed haploid; missing haplotig = homozygous |
| 3 | 02 | diagnostic definition | allele private to one genotype's allele set | Ours | heterozygous sites usable when a private allele exists |
| 4 | 04/11 | alignment | STAR two-pass, uniquely mapped reads only | Ours (standard practice) | avoid multi-mapping paralog reads |
| 5 | 06/08/12 | counting | bcftools mpileup `-q20 -Q20 -B -d10000`, no genotype calling | Ours; same approach (raw AD counts into baymobil) as Paajanen et al. 2025 | calling would discard low-fraction donor reads |
| 6 | 07 | positions counted | genome exonic markers ∪ RNA positions passing L1 in ≥2 samples | Ours | keeps counting tractable |
| 7 | 09 L1 | per-sample allele support | depth ≥10, allele reads ≥3, fraction ≥0.10, both orientations if ≥4 reads | Ours | discovery standard for homografts |
| 8 | 09 L2 | reproducibility | ≥2 replicates per genotype × tissue | Ours | removes one-off library errors |
| 9 | 09 L3 | RNA editing | drop A→G (+ genes) / T→C (− genes) without genome support | Ours | A-to-I editing signature |
| 10 | 09 L4 | donor-private (RNA-only markers) | donor allele in any donor tissue; recipient depth ≥10 in ≥2 reps; recipient fraction ≤2% | Ours | v3: any donor tissue (in SR grafts the Ruggeri part is the root) |
| 11 | 13 | marker tested | ≥1 heterograft read | Ours | test every expressed marker, no homograft-expression requirement |
| 12 | 13 | recipient-homograft removal | donor allele in recipient homograft ≥1 read and >2% → remove | Ours | implements "− recipient self-graft SNPs" |
| 13 | 13 | variant type | SNPs only | Ours | indels add ~2% testable genes; higher RNA-seq error; not covered by baymobil paper |
| 14 | 13 | Bayes factor | baymobil exact Beta-binomial, nmax = 10 | **Published** (Tomkins et al. 2022; baymobil) | nmax 10 = package default |
| 15 | 13 | error-rate data | Nh1 recipient homograft same tissue; Nh2 donor homograft all tissues | Published structure (two homografts); tissue choices Ours | donor all tissues gives coverage for genes not expressed in the donor tissue of the same name |
| 16 | 13 | homograft depth cap | 2000 (counts scaled) | Ours | computation time; note: stronger cap → stronger n=0 penalty (≈ −3.1 vs −2.5 at 300) |
| 17 | 13b | gene score | sum of per-SNP log10BF | **Published** (Tomkins et al. 2022) | penalises genes with many SNPs and low mobile fraction (conservative) |
| 18 | 13b | evidence threshold | gene log10BF ≥ 1 | **Published** (Tomkins 2022; Paajanen 2025 use log BF > 1) | |
| 19 | 13b | replicate rule | ≥2 replicates pass (main); ≥3 also reported | Ours | concordance check: 396 observed vs 28.8 expected by chance |
| 20 | 13b | gene-pooled BF | counts pooled over SNPs ≥300 bp apart, one BF | Ours (not published) | sensitivity analysis only |
| 21 | all | excluded sample | SRR4 | Ours | labelled Ruggeri root but mostly Shiraz RNA (label error / mix-up) |
| 22 | 18 | high-confidence subset | obs/exp ≥5 vs error + expression-proportional background, both sources | Ours (arbitrary) | prioritisation for wet-lab validation, NOT a correction of the main result |
