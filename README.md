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

**Step 2 — RNA confirmation and homograft filtering.**
Align the two self-graft (homograft) RNA-seq sets (Shz/Shz, Ru/Ru) to `T2T_ref` with
**STAR** [5]. At each diagnostic position, count ref/alt reads with
**GATK ASEReadCounter** [6] (or **BCFtools** mpileup [3]). Keep positions that are
(a) expressed, and (b) show **only the self allele** in the matching homograft; discard
positions that leak the other allele (mis-mapping / paralogues, or a heterozygous site
whose haplotig was missed in the assembly). Optionally score each site with a
beta-binomial model in **R** [7] / **tidyverse** [8].

**Step 3 — Master diagnostic database.**
Validated, expressed, clean diagnostic positions, grouped per gene, prioritising genes
with >=2 diagnostic sites within one amplicon and preferring indels (indels give a size
difference usable by gel without a restriction enzyme).

**Step 4 — Mobile-transcript detection in heterografts.**
Align heterograft RNA-seq (e.g. Shz-scion/Ru-root) to `T2T_ref` with STAR [5], and at
each diagnostic position count reads carrying the **donor allele** (whatever it is, ref
or alt) with ASEReadCounter [6]. Require **multiple diagnostic sites per transcript to
agree** (co-occurring evidence) and compare to the homograft background.

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
├── docs/                    <- experiment_background.md, analysis_log.md, method_notes.md
├── metadata/                <- sample_sheet.tsv, mobile_candidates.tsv
├── scripts/                 <- 00_setup_env.sh, 01_dipcall.sh, 02_build_diagnostic_sites.R, tests/
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
| GATK (ASEReadCounter) | allele-specific read counts | McKenna et al. 2010, *Genome Res* 20:1297–1303. doi:10.1101/gr.107524.110 [6] |
| R | statistics | R Core Team (version used). https://www.R-project.org/ [7] |
| tidyverse / data.table | data handling in R | Wickham et al. 2019, *JOSS* 4:1686. doi:10.21105/joss.01686 [8] |
| AdapterRemoval v2 | read trimming (upstream) | Schubert et al. 2016, *BMC Res Notes* 9:88. doi:10.1186/s13104-016-1900-2 [9] |
| gffread | GFF/GTF handling, transcript extraction | Pertea & Pertea 2020, *F1000Res* 9:304. doi:10.12688/f1000research.23297.2 [10] |
| SeqKit | FASTA/FASTQ manipulation | Shen et al. 2016, *PLoS ONE* 11:e0163962. doi:10.1371/journal.pone.0163962 [11] |
| featureCounts (Subread) | gene-level counts (expression) | Liao et al. 2014, *Bioinformatics* 30:923–930. doi:10.1093/bioinformatics/btt656 [12] |
| FastQC / MultiQC | QC reports | Ewels et al. 2016, *Bioinformatics* 32:3047–3048. doi:10.1093/bioinformatics/btw354 [13] |

> Record the EXACT version of every tool you run (from `environment.yml` /
> `conda list`) in `docs/analysis_log.md`, and cite these versions in any thesis/paper.

## 7. Links
- Thesis repo (frozen): `AdelaideUniPlantMobileRNA/Mobile_RNA_Grapevine_Mapping_Based` @ `thesis-submission-2026-10`
- Contact: Fei Zheng — zhengfei198905@gmail.com
