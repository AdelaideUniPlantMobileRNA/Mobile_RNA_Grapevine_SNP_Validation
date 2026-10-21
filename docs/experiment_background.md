# Experiment background (biological + data context)

## 1. Biological question & hypothesis
- **Question:** Which transcripts move across the Shiraz / 140-Ruggeri graft junction,
  and can they be reliably identified and validated using genotype-distinguishing SNPs?
- **Hypothesis:** A subset of transcripts is graft-mobile between rootstock and scion;
  these can be detected by diagnostic SNP/indel differences between the two genotypes.

## 2. Plant material
- **Scion / rootstock genotypes:** Shiraz (*Vitis vinifera* L.) and 140-Ruggeri
  (*V. berlandieri* x *V. rupestris*).
- **Graft combinations (4):**
  - Shiraz scion / Ruggeri rootstock  (heterograft)
  - Ruggeri scion / Shiraz rootstock  (reciprocal heterograft)
  - Shiraz scion / Shiraz rootstock   (homograft control)
  - Ruggeri scion / Ruggeri rootstock (homograft control)

## 3. Sample layout (RNA-seq)
- **Tissues (4):** leaf, bud (latent compound bud), phloem (young stem scraping), root.
- **Replicates:** 4 per sample, EXCEPT Ruggeri/Ruggeri homograft = 3 replicates.
- **Total:** ~60 RNA-seq libraries.
- Illumina paired-end; reads already adapter/quality trimmed (AdapterRemoval v2).

### Data location (HPC, Phoenix)
- Trimmed FASTQ: `/scratchdata1/users/a1709557/Own_RNA/02.TrimmedData/T1T2`

## 4. Reference genomes and annotations
### Genome assemblies (FASTA)
| Genotype | Haplotype | File | Source |
|---|---|---|---|
| PN40024 (V. vinifera) T2T | ref | `T2T_ref.fasta` | grapedia.org |
| PN40024 T2T | alt | `T2T_alt.fasta` | grapedia.org |
| Shiraz | primary | `SHZ_1654_primary.fasta` | doi:10.1155/2023/6686706 |
| Shiraz | haplotigs | `SHZ_1654_haplotigs.fasta` | doi:10.1155/2023/6686706 |
| 140-Ruggeri | primary | `Primary_contigs_140Ruggeri.fasta` | our lab (Zhou-Tsang et al. 2022, unpubl., with permission) |
| 140-Ruggeri | alt | `alt_contigs_140Ruggeri.fasta` | our lab |

### Annotations
| File | Applies to | Format |
|---|---|---|
| `PN40024_5.1_on_T2T_ref_with_names.gff3` | PN40024 T2T | GFF3 (grapedia.org) |
| `SHZ_1654_primary_annotation.txt` | Shiraz primary | txt |
| `SHZ_1654_haplotigs_annotation.txt` | Shiraz haplotigs | txt |
| `140R.P_contigs.gtf` | Ruggeri primary | GTF |
| `140R.alt_contigs.gtf` | Ruggeri alt | GTF |

**Note on "two FASTA + one gene ID on two contigs":**
Each genotype is a **phased diploid assembly**: the two FASTA files are its two
**haplotypes**. The same gene ID appearing on two different contigs (e.g.
`SHZ1654_000001` on a primary contig and on a haplotig) is the **two allelic copies**
of one gene (one per haplotype). This is expected. For diagnostic-SNP calling, both
haplotypes are used together (dipcall treats them as the diploid). For gene-level
expression, `featureCounts -g gene_id` collapses the two allelic copies into one gene.

## 5. Known issues / caveats
- 140-Ruggeri homograft has only **3** replicates (others have 4).
- 140-Ruggeri is **not V. vinifera**; alignment to the PN40024 T2T backbone leaves
  blind spots in divergent regions.
- Phloem samples are **partial tissue** (stem scraping), i.e. a snapshot of transiting RNA.
- Both genotypes are **highly heterozygous clones**; only homozygous-but-different sites
  are usable as clean diagnostic markers.

## 6. Key references
- Shiraz genome: doi:10.1155/2023/6686706
- PN40024 T2T: Shi et al. 2023 (grapedia.org); PN40024.5.1 annotation.
- 140-Ruggeri genome: Zhou-Tsang et al. 2022 (poster; unpublished; used with permission).
- Prior method (thesis): Mobile_RNA_Grapevine_Mapping_Based.
