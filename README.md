# Mobile_RNA_Grapevine_SNP_Validation

SNP-based identification and validation of graft-mobile RNA in Shiraz / 140-Ruggeri
grapevine grafts. Follow-up to **Mobile_RNA_Grapevine_Mapping_Based**
(PhD thesis Chapter 4; frozen at tag `thesis-submission-2026-10`).

---

## 1. Goal
Identify graft-mobile transcripts using a **SNP-based**, genotype-resolved method,
and build a **validated diagnostic SNP/indel database** that distinguishes Shiraz
from 140-Ruggeri, suitable for laboratory validation (CAPS RT-PCR; Nanopore long-read).

## 2. Rationale
Both Shiraz and 140-Ruggeri are **highly heterozygous**, so naive SNP calling cannot
cleanly assign a read to a genotype. The key design choice is to **co-genotype both
genotypes on a single common coordinate backbone (PN40024 T2T)** and then keep only
positions carrying a **private allele** (an allele present in one genotype and absent
from the other). A read carrying a donor-private allele in the recipient tissue is
evidence of mobility.

> Note: 140-Ruggeri is an interspecific hybrid (*V. berlandieri* x *V. rupestris*),
> not *V. vinifera*. Aligning it to the *V. vinifera* PN40024 backbone leaves blind
> spots in Ruggeri-divergent regions (see Limitations).

## 3. Workflow

**Step 1 — Genome-based diagnostic catalogue (co-genotyping on T2T).** *(implemented — see `scripts/` and `docs/method_notes.md`)*
Each genotype is a diploid represented by two haplotype assemblies (primary +
haplotig/alt). Align both haplotypes of each genotype to `T2T_ref` and call a diploid
VCF + per-haplotype confident BED with **dipcall** [1] (which uses **minimap2** [2];
Ruggeri is relaxed to the `asm20` preset for its greater divergence). Normalise and
merge the two VCFs on T2T coordinates with **BCFtools** [3].

Because both assemblies were processed with Purge Haplotigs allelic-contig
reassignment, the **primary assembly is a collapsed haploid representation** and the
haplotig/alt contigs cover only the heterozygous fraction. The callable region used
is therefore each genotype's **primary (hap1) coverage** — *not* dipcall's `dip.bed`,
which requires both haplotypes to align and discards the homozygous majority. Within
primary coverage, a missing second haplotype is read as **homozygous**, and a no-call
is read as confident **0/0**.

Each genotype is reduced to an **allele set**, and a position is **diagnostic** when
one genotype carries an allele the other lacks (**directional private alleles**):
`shiraz_private` marks scion→rootstock mobility, `ruggeri_private` marks
rootstock→scion mobility. Both homozygous-but-different sites *and* heterozygous sites
contribute, as long as a private allele exists (interval operations with **BEDTools**
[4]). Sites are tagged with gene/exon context from the PN40024-on-T2T GFF3.

**Step 2 — RNA evidence from homografts (marker discovery + confirmation).** *(implemented — `scripts/03–09`, `docs/step2_notes.md`)*
Homograft RNA-seq (SS = Shiraz/Shiraz, RR = 140-Ru/140-Ru; 4 tissues: B, L, PY, R) is
aligned to `T2T_ref` with **STAR** [5] (two-pass, uniquely mapped reads only).
Alternate alleles are discovered in exons and allele-specific read counts
(AD/ADF/ADR) are taken at all candidate positions with **BCFtools mpileup** [3]
(no genotype calling; GATK ASEReadCounter was tested and dropped because it rejects
positions with more than one variant context). Hard-threshold filters L1–L4
(`docs/assumptions.md`) define an **RNA-only** marker set that adds donor alleles
missed by the assemblies (assembly gaps, clone-specific variants); a donor allele is
accepted if it is seen in **any** donor tissue and is absent (≤2%) from the
tissue-matched recipient homograft.

**Step 3 — Marker set used for detection.**
Genome markers (Step 1: exonic, donor-private SNP alleles; used directly, no RNA
support required) ∪ RNA-only markers (Step 2), minus any marker whose donor allele is
seen in the tissue-matched **recipient** homograft (≥1 read and >2%). SNPs only;
indels were evaluated and not used (they would add ~2% testable genes; see
`docs/analysis_log.md`).

**Step 4 — Mobile-transcript detection in heterografts (Bayesian).** *(implemented — `scripts/11–15`)*
Heterograft RNA-seq is aligned the same way; donor- and total-allele read counts are
taken at every marker. Every marker with ≥1 heterograft read is tested with
**baymobil** [14] (exact Bayesian inference, Tomkins et al. 2022 [15]): the error rate
is estimated from the tissue-matched recipient homograft (local) and from the donor
homograft (all tissues pooled), and a per-SNP log10 Bayes factor compares "all donor
reads are errors" with "some reads are mobile". Per-SNP log10BF are **summed per gene**
(published scheme [15]); a replicate passes if the gene log10BF ≥ 1, and a
gene × tissue is called mobile if it passes in **≥2 biological replicates**. A
gene-pooled Bayes factor (our extension, not published) is reported side by side.

**Step 5 — Verification.** *(scripts/16–19, check_gene_labels.py)*
Global donor-allele fraction heterograft vs homograft; replicate concordance vs
chance; heterograft-specific contamination check following Paajanen et al. 2025 [16];
independent re-counting with samtools and recomputation with baymobil's `run_bayes`
on randomly chosen markers; gene-ID-vs-coordinate check against the GFF3.

## 4. Known limitations
- **Ruggeri divergence:** Ruggeri-specific / high-divergence regions do not anchor to
  the *V. vinifera* T2T backbone and are blind spots (for both catalogue building and
  read recovery).
- **Heterozygous sites carry less information:** at a heterozygous site only the
  allele that is private to one genotype is diagnostic; the shared allele is
  uninformative. Such sites are still used when a private allele exists, but give
  one-directional evidence only.
- **Hidden-heterozygosity assumption:** treating an absent second haplotype as
  homozygous assumes haplotigs captured all heterozygous sites; an assembly gap could
  mis-read a true heterozygote as homozygous. The Step-2 homograft RNA filter catches
  these empirically (see `docs/method_notes.md`).
- **Irreducible blind spots:** where donor and recipient carry the same allele along a
  whole transcript, origin cannot be determined by any SNP method.
- **Absence of evidence:** low-abundance mobile transcripts may be below detection depth
  (a negative is "not detected", not "absent").

## 5. Directory structure
```
.
├── README.md
├── environment.yml          <- conda env (command-line tools)
├── config/config.yaml       <- all paths & parameters
├── docs/                    <- experiment_background.md, analysis_log.md, method_notes.md,
│                               step2_notes.md (Steps 2–5), assumptions.md
├── metadata/                <- sample_sheet.tsv, mobile_candidates.tsv
├── scripts/                 <- 00–19 pipeline scripts (run order: docs/step2_notes.md), tests/
├── data/ (not in git)       <- raw/ processed/ external/
├── results/ (not in git)
└── logs/ (not in git)
```
See `docs/method_notes.md` for the full reasoning behind the allele-handling choices
(primary-coverage callable region, the Purge-Haplotigs convention, how each region
type is treated).

## 6. Software and citations
| Tool | Use | Citation |
|---|---|---|
| dipcall | diploid assembly → VCF vs T2T | Li et al. 2018, *Nat Methods* 15:595–597. doi:10.1038/s41592-018-0054-7 [1] |
| minimap2 | assembly/long-read alignment | Li 2018, *Bioinformatics* 34:3094–3100. doi:10.1093/bioinformatics/bty191 [2] |
| SAMtools / BCFtools | VCF norm/merge, pileup, BAM ops | Danecek et al. 2021, *GigaScience* 10:giab008. doi:10.1093/gigascience/giab008 [3] |
| BEDTools | interval operations | Quinlan & Hall 2010, *Bioinformatics* 26:841–842. doi:10.1093/bioinformatics/btq033 [4] |
| STAR | RNA-seq alignment | Dobin et al. 2013, *Bioinformatics* 29:15–21. doi:10.1093/bioinformatics/bts635 [5] |
| GATK (ASEReadCounter) | tested, NOT used (see Step 2) | McKenna et al. 2010, *Genome Res* 20:1297–1303. doi:10.1101/gr.107524.110 [6] |
| R | statistics | R Core Team (version used). https://www.R-project.org/ [7] |
| tidyverse / data.table | data handling in R | Wickham et al. 2019, *JOSS* 4:1686. doi:10.21105/joss.01686 [8] |
| AdapterRemoval v2 | read trimming (upstream) | Schubert et al. 2016, *BMC Res Notes* 9:88. doi:10.1186/s13104-016-1900-2 [9] |
| gffread | GFF/GTF handling, transcript extraction | Pertea & Pertea 2020, *F1000Res* 9:304. doi:10.12688/f1000research.23297.2 [10] |
| SeqKit | FASTA/FASTQ manipulation | Shen et al. 2016, *PLoS ONE* 11:e0163962. doi:10.1371/journal.pone.0163962 [11] |
| featureCounts (Subread) | gene-level counts (expression) | Liao et al. 2014, *Bioinformatics* 30:923–930. doi:10.1093/bioinformatics/btt656 [12] |
| FastQC / MultiQC | QC reports | Ewels et al. 2016, *Bioinformatics* 32:3047–3048. doi:10.1093/bioinformatics/btw354 [13] |
| baymobil | Bayesian mobility test | Hoerbst, Morris & Tomkins, Research Square preprint 2023. doi:10.21203/rs.3.rs-2520491/v1 [14] |
| (method) | exact Bayesian inference for graft-mobile transcripts | Tomkins et al. 2022, *J R Soc Interface* 19:20220644. doi:10.1098/rsif.2022.0644 [15] |
| (method) | re-analysis / recommended checks for mobile mRNA | Paajanen et al. 2025, *Nat Plants* 11:977–984. doi:10.1038/s41477-025-01979-x [16] |
| Python: pandas / NumPy | data handling in Python | McKinney 2010 (SciPy proc.); Harris et al. 2020, *Nature* 585:357–362. doi:10.1038/s41586-020-2649-2 [17] |

> Record the EXACT version of every tool you run (from `environment.yml` /
> `conda list`) in `docs/analysis_log.md`, and cite these versions in any thesis/paper.

## 7. Links
- Thesis repo (frozen): `AdelaideUniPlantMobileRNA/Mobile_RNA_Grapevine_Mapping_Based` @ `thesis-submission-2026-10`
- Contact: Fei Zheng — zhengfei198905@gmail.com
