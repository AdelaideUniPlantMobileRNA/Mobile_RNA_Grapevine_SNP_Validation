#!/usr/bin/env Rscript
# =====================================================================
# Step 1 (cont.): build the CANDIDATE genotype-distinguishing SNP/indel
# database between Shiraz and Ruggeri, from the merged dipcall genotypes.
#
# KEY ASSEMBLY FACT (Onetto et al. 2023, doi:10.1155/2023/6686706):
#   Both assemblies were processed with Purge Haplotigs ("allelic contig
#   reassignment"). The PRIMARY assembly is a collapsed haploid
#   representation; the alternate copy (haplotig / alt-contig) is only
#   kept where the two haplotypes DIFFER (heterozygous/divergent regions).
#   => Where the second haplotype is ABSENT, the genotype is HOMOZYGOUS
#      (or hemizygous) for the primary allele -- NOT "unknown".
#
# Consequences for this script:
#   * Callable region = PRIMARY (hap1) coverage BED  (shiraz.hap1.bed /
#     ruggeri.hap1.bed), NOT dipcall's dip.bed. dip.bed requires BOTH
#     haplotypes to align and is far too conservative for a purged
#     primary assembly (it would discard the homozygous majority).
#   * A half-call "X|." (hap1 called, hap2 missing) is read as HOMOZYGOUS
#     for allele X. In allele-set terms it is simply a one-element set.
#   * A no-call "./." inside the primary (hap1) coverage = 0/0 (confident
#     reference). A no-call OUTSIDE coverage = unknown => position dropped.
#
# Diagnostic logic = ALLELE SETS + DIRECTION (private alleles):
#   shiraz_private  = alleles Shiraz carries that Ruggeri does NOT
#                     -> markers to detect SCION (Shiraz) RNA in ROOTSTOCK
#   ruggeri_private = alleles Ruggeri carries that Shiraz does NOT
#                     -> markers to detect ROOTSTOCK (Ruggeri) RNA in SCION
#   A site is a candidate if EITHER private set is non-empty.
#
# These are CANDIDATES. Step 2 (homograft SS/RR RNA allele counts)
# empirically confirms each genotype's real expressed allele set and
# filters mis-mapping, especially catching any heterozygous site whose
# haplotig was missed (which would otherwise look homozygous here).
#
# Usage:
#   Rscript 02_build_diagnostic_sites.R <outdir> [gff3]
#     <outdir> : dir holding merged_genotypes.tsv + *.hap1.bed   (default results/dipcall)
#     [gff3]   : optional gene annotation on the T2T backbone for gene/exon tagging
# =====================================================================
suppressPackageStartupMessages(library(data.table))
setDTthreads(0)

args    <- commandArgs(trailingOnly = TRUE)
outdir  <- if (length(args) >= 1) args[1] else "results/dipcall"
gffpath <- if (length(args) >= 2) args[2] else NA_character_

gt <- fread(file.path(outdir, "merged_genotypes.tsv"), header = FALSE,
            sep = "\t", quote = "",
            col.names = c("CHROM","POS","REF","ALT","shiraz_GT","ruggeri_GT"))
gt[, CHROM := as.character(CHROM)][, POS := as.integer(POS)]
for (cc in c("REF","ALT","shiraz_GT","ruggeri_GT"))
  set(gt, j = cc, value = as.character(gt[[cc]]))
message("Merged positions read: ", nrow(gt))

# ---- primary (hap1) coverage: used ONLY to infer 0/0 where a genotype made no call ----
read_bed <- function(f) {
  b <- fread(f, header = FALSE, select = 1:3, col.names = c("chrom","start","end"))
  b[, start := start + 1L]            # BED 0-based half-open -> 1-based inclusive
  setkey(b, chrom, start, end)
  b[]
}
shz_bed <- read_bed(file.path(outdir, "shiraz.hap1.bed"))
rug_bed <- read_bed(file.path(outdir, "ruggeri.hap1.bed"))

q <- data.table(chrom = gt$CHROM, start = gt$POS, end = gt$POS)
setkey(q, chrom, start, end)
in_shz <- !is.na(foverlaps(q, shz_bed, type = "any", mult = "first", which = TRUE))
in_rug <- !is.na(foverlaps(q, rug_bed, type = "any", mult = "first", which = TRUE))

# ---- split phased GT into two haplotype tokens ("1|." -> "1",".") ----
gt[, c("s1","s2") := tstrsplit(shiraz_GT,  "[|/]")]
gt[, c("r1","r2") := tstrsplit(ruggeri_GT, "[|/]")]

# ---- map a GT index token to its actual allele base ----
first_alt <- sub(",.*", "", gt$ALT)                 # fast path: ALT index 1
altlist   <- strsplit(gt$ALT, ",", fixed = TRUE)     # slow path: multiallelic only
tok2base <- function(tok) {
  out <- rep(NA_character_, length(tok))
  i0  <- which(tok == "0");  out[i0] <- gt$REF[i0]   # 0 -> REF
  i1  <- which(tok == "1");  out[i1] <- first_alt[i1] # 1 -> first ALT
  ihi <- which(!is.na(tok) & !(tok %in% c("0","1",".")))  # >=2 -> nth ALT (rare)
  if (length(ihi)) {
    k <- as.integer(tok[ihi])
    out[ihi] <- vapply(seq_along(ihi), function(j) altlist[[ihi[j]]][k[j]], character(1))
  }
  out
}
gt[, sa1 := tok2base(s1)][, sa2 := tok2base(s2)]
gt[, ra1 := tok2base(r1)][, ra2 := tok2base(r2)]

# ---- infer confident homozygous reference (0/0) for no-call inside primary coverage ----
gt[is.na(sa1) & is.na(sa2) & in_shz, sa1 := REF]
gt[is.na(ra1) & is.na(ra2) & in_rug, ra1 := REF]

# ---- drop positions still unknown for EITHER genotype (outside its primary coverage) ----
gt <- gt[!(is.na(sa1) & is.na(sa2)) & !(is.na(ra1) & is.na(ra2))]
message("Positions callable in BOTH genotypes: ", nrow(gt))

# ---- canonicalise each genotype to an ordered allele pair A1[,A2] ----
# Handles any slot order (incl. ".|1") and collapses X|X and X|. to one allele:
#   A1 = the (always-present) first allele;  A2 = the distinct second allele, else NA.
canon1 <- function(a1, a2) fcoalesce(a1, a2)
canon2 <- function(a1, a2) fifelse(!is.na(a1) & !is.na(a2) & a1 != a2, a2, NA_character_)
gt[, `:=`(sA1 = canon1(sa1, sa2), sA2 = canon2(sa1, sa2),
          rA1 = canon1(ra1, ra2), rA2 = canon2(ra1, ra2))]

# ---- directional private alleles (each side has at most 2 distinct alleles) ----
inset <- function(x, b1, b2)                         # is base x present in the set {b1,b2}?
  (!is.na(x)) & (((x == b1) & !is.na(b1)) | ((x == b2) & !is.na(b2)))

gt[, sA1_priv := !is.na(sA1) & !inset(sA1, rA1, rA2)]
gt[, sA2_priv := !is.na(sA2) & !inset(sA2, rA1, rA2)]
gt[, rA1_priv := !is.na(rA1) & !inset(rA1, sA1, sA2)]
gt[, rA2_priv := !is.na(rA2) & !inset(rA2, sA1, sA2)]

join_two <- function(a1, keep1, a2, keep2) {         # comma-join the kept alleles
  v1 <- fifelse(keep1, a1, NA_character_)
  v2 <- fifelse(keep2, a2, NA_character_)
  fifelse(!is.na(v1) & !is.na(v2), paste(v1, v2, sep = ","),
   fifelse(!is.na(v1), v1, fifelse(!is.na(v2), v2, NA_character_)))
}
allele_set <- function(a1, a2) fifelse(is.na(a2), a1, paste(a1, a2, sep = ","))

gt[, shiraz_alleles   := allele_set(sA1, sA2)]
gt[, ruggeri_alleles  := allele_set(rA1, rA2)]
gt[, shiraz_private   := join_two(sA1, sA1_priv, sA2, sA2_priv)]
gt[, ruggeri_private  := join_two(rA1, rA1_priv, rA2, rA2_priv)]

# ---- keep candidates: at least one side carries a private allele ----
diag <- gt[!is.na(shiraz_private) | !is.na(ruggeri_private)]
diag[, detects_shiraz_in_rootstock := !is.na(shiraz_private)]  # scion -> rootstock marker
diag[, detects_ruggeri_in_scion    := !is.na(ruggeri_private)] # rootstock -> scion marker
diag[, is_indel := (nchar(REF) != nchar(sA1)) |
                   (!is.na(sA2) & nchar(REF) != nchar(sA2)) |
                   (nchar(REF) != nchar(rA1)) |
                   (!is.na(rA2) & nchar(REF) != nchar(rA2))]

# ---- optional: tag gene / exon from a T2T-backbone GFF3 ----
if (!is.na(gffpath) && file.exists(gffpath)) {
  g <- fread(cmd = paste("grep -v '^#'", shQuote(gffpath)),
             header = FALSE, sep = "\t", quote = "")
  setnames(g, 1:9, c("chrom","src","type","start","end","score","strand","frame","attr"))
  g[, chrom := as.character(chrom)]
  genes <- g[type == "gene"]
  genes[, gene_id   := sub(".*ID=([^;]+).*",   "\\1", attr)]
  genes[, gene_name := fifelse(grepl("Name=", attr),
                               sub(".*Name=([^;]+).*", "\\1", attr), NA_character_)]
  genes <- genes[, .(chrom, start, end, gene_id, gene_name)]
  exons <- g[type %in% c("exon","CDS"), .(chrom, start, end)]
  setkey(genes, chrom, start, end); setkey(exons, chrom, start, end)

  # sanity: do GFF and VCF use the same chromosome names?
  message("GFF: ", nrow(genes), " genes, ", nrow(exons), " exon/CDS records.")
  message("  GFF chroms (head): ", paste(head(unique(genes$chrom)), collapse = ","))
  message("  VCF chroms (head): ", paste(head(unique(diag$CHROM)), collapse = ","))
  shared <- intersect(unique(genes$chrom), unique(diag$CHROM))
  message("  shared chrom names: ", length(shared),
          if (length(shared) == 0) "  <-- NAME MISMATCH: gene/exon will be empty" else "")

  dq <- data.table(chrom = diag$CHROM, start = diag$POS, end = diag$POS)
  setkey(dq, chrom, start, end)
  gi <- foverlaps(dq, genes, type = "any", mult = "first")
  diag[, gene_id   := gi$gene_id]
  diag[, gene_name := gi$gene_name]
  diag[, is_exonic := !is.na(foverlaps(dq, exons, type = "any", mult = "first", which = TRUE))]
} else {
  diag[, `:=`(gene_id = NA_character_, gene_name = NA_character_, is_exonic = NA)]
  if (is.na(gffpath)) message("No GFF given -> gene/exon columns left NA.")
}

out <- diag[, .(CHROM, POS, REF, ALT,
                shiraz_alleles, ruggeri_alleles,
                shiraz_private, ruggeri_private,
                detects_shiraz_in_rootstock, detects_ruggeri_in_scion,
                is_indel, gene_id, gene_name, is_exonic)]
setorder(out, CHROM, POS)
fwrite(out, file.path(outdir, "diagnostic_sites.tsv"),
       sep = "\t", na = "NA", quote = FALSE)

message("\nCandidate diagnostic sites: ", nrow(out),
        "  (SNPs: ", sum(!out$is_indel), ", indels: ", sum(out$is_indel), ")")
message("  detect Shiraz->rootstock (shiraz-private): ", sum(out$detects_shiraz_in_rootstock))
message("  detect Ruggeri->scion   (ruggeri-private): ", sum(out$detects_ruggeri_in_scion))
if (!all(is.na(out$is_exonic)))
  message("  in a gene: ", sum(!is.na(out$gene_id)),
          "   exonic: ", sum(out$is_exonic, na.rm = TRUE),
          "   (genic+exonic are the Step-2 priority set)")
message("  -> ", file.path(outdir, "diagnostic_sites.tsv"))
