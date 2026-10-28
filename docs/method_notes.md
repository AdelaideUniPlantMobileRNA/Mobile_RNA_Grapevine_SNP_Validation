# Method notes — how diagnostic sites are defined

This explains the allele-handling choices in Step 1.2
(`02_build_diagnostic_sites.R`). It is the reasoning record for the pipeline
and for the thesis methods.

## Starting point

Positions enter `merged_genotypes.tsv` only when **at least one** genotype
differs from the T2T reference. Positions where Shiraz and Ruggeri are both
identical to T2T never appear — they carry no information and are correctly
absent. So every candidate position already has some difference somewhere; the
job of Step 1.2 is to decide, per position, whether Shiraz and Ruggeri differ
**from each other**, and in which direction.

## Callable region = primary (hap1) coverage, not dipcall's `dip.bed`

Both assemblies were processed with **Purge Haplotigs** ("allelic contig
reassignment"; Onetto et al. 2023, doi:10.1155/2023/6686706). As the Shiraz
assembler describes it: the assembly is split into two contigs **only where the
two haplotypes differ** (heterozygous/divergent regions); where they are the
same (homozygous), a single primary contig is kept. The **primary** assembly is
therefore a collapsed haploid representation of the whole genome, and the
**haplotig / alt** contigs cover only the heterozygous fraction.

dipcall's `dip.bed` marks regions where **both** haplotypes aligned. For a
purged primary assembly that discards the homozygous majority of the genome
(for Shiraz it was ~36 Mb out of a ~187 Mb primary). That is wrong for this
data: where the second haplotype is absent, the genotype is **homozygous**, not
"unknown". So we use each genotype's **primary (hap1) coverage** as the callable
region (Shiraz ~187 Mb, Ruggeri ~213 Mb).

Consistency check: Ruggeri is an interspecific hybrid (*V. berlandieri* ×
*V. rupestris*) and so is heterozygous across almost the whole genome — its
haplotigs are nearly as large as its primary (hap2 ≈ hap1). Shiraz, a *vinifera*
cultivar, is less heterozygous, so its haplotigs cover far less. The smaller
Shiraz haplotig coverage is therefore *expected homozygosity*, not a poor
assembly — exactly matching the assembler's account.

## Genotype interpretation

Each genotype's two phased haplotype calls are reduced to an **allele set**:

| dipcall genotype | meaning here | allele set |
|------------------|--------------|------------|
| `0\|1`, `1\|0` | heterozygous (both haplotypes aligned) | {REF, ALT} |
| `1\|1`, `0\|0` | homozygous | {that allele} |
| `1\|.`, `.\|1`, `0\|.` | second haplotype absent → **homozygous** (primary-assembly convention) | {primary allele} |
| `./.` **inside** primary coverage | no variant, genome aligned → reference | {REF} |
| `./.` **outside** primary coverage | not callable | unknown → position dropped |
| multiallelic (e.g. `1\|1` vs `2\|2`, ALT=`A,T`) | index mapped to the correct ALT base | {correct base} |

A one-element allele set *is* homozygosity, so `X|.` needs no special case — it
simply yields `{X}`.

## What counts as diagnostic (directional private alleles)

For each position we compute, per direction, the alleles one genotype carries
that the other lacks:

* `shiraz_private  = shiraz_alleles \ ruggeri_alleles`
  → a marker to detect **scion (Shiraz) RNA appearing in rootstock tissue**.
* `ruggeri_private = ruggeri_alleles \ shiraz_alleles`
  → a marker to detect **rootstock (Ruggeri) RNA appearing in scion tissue**.

A position is kept if **either** private set is non-empty.

### The three region types asked about

1. **Highly heterozygous regions** (both haplotypes assembled, e.g. `0|1`).
   Included. Example: T2T=A, Shiraz `0|1`={A,G}, Ruggeri `1|1`={G}.
   `shiraz_private = {A,G}\{G} = G` (usable in the Shiraz→rootstock direction);
   `ruggeri_private = {G}\{A,G} = ∅`. Heterozygous sites carry less information
   (the shared allele is uninformative) but are still valid whenever one
   genotype has an allele the other lacks.

2. **Shiraz differs from T2T but Ruggeri equals T2T.** Included, not lost.
   Example: T2T=A, Shiraz `1|1`={C}, Ruggeri has no variant record. The position
   is in the file because of Shiraz's variant; Ruggeri's missing call lies in
   its primary coverage and is filled as `0/0`={A}. Result:
   `shiraz_private=C`, `ruggeri_private=A` — diagnostic both ways.

3. **Both homozygous (the Purge-Haplotigs convention).** Handled correctly:
   * homozygous for the **same** allele → no private allele → dropped (nothing
     to distinguish — correct);
   * homozygous for **different** alleles (e.g. Shiraz {T}, Ruggeri {C}) → both
     sides private → diagnostic, and the cleanest marker class (both
     homozygous, disjoint alleles — easiest to design RT-PCR against).

## The one residual assumption, and how Step 2 covers it

Treating "haplotig absent" as homozygous assumes the haplotigs captured **all**
truly heterozygous sites. If an assembly gap left a genuinely heterozygous
region with only the primary aligned, it would show as `X|.` and be mis-read as
homozygous. The directional risk is asymmetric:

* mis-calling a **Shiraz** het as homozygous only **loses** a marker
  (false negative — the primary allele it reports is still genuinely Shiraz's);
* mis-calling a **Ruggeri** (recipient) het as homozygous could create a
  **false-positive** "Shiraz-private" call. Ruggeri's haplotigs are near
  complete, so this risk is low.

**Step 2 (homograft RNA) is the empirical safeguard.** Self-graft Shiraz (SS)
and Ruggeri (RR) RNA give each genotype's real expressed allele set: if a
candidate site shows two alleles in the self-graft, it was actually
heterozygous and the corresponding marker is discarded. The same RNA step also
removes sites distorted by read mis-mapping. Only sites that survive this
confirmation are used for mobility calls.
