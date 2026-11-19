#!/usr/bin/env Rscript
# =====================================================================
# Step 2.6 (v3): build the master diagnostic marker database from homograft RNA.
#
# v3 change (only L4): a donor allele is accepted if it passes L1-L3 in ANY tissue
#   of the donor self-graft (e.g. RR root, RR leaf, ...), not only in the same tissue
#   as the recipient. It is then checked against the recipient self-graft of EACH
#   tissue (tissue-matched background, as before). Reason: in an SR graft the Ruggeri
#   part is the root; a Ruggeri root mRNA moving to Shiraz leaf need not be expressed
#   in Ruggeri leaf. Output: master_diagnostic_db_v3.tsv (v1 file is not touched).
#
# Two DIFFERENT standards apply in this project:
#   * DISCOVERY (here): STRICT — define, per genotype, which alleles are really
#     expressed, excluding sequencing error, contamination, strand artefacts and
#     RNA editing. This is the trusted "truth" database.
#   * DETECTION (Step 4, heterograft): SENSITIVE — look for a few donor-allele
#     reads turning up in the recipient tissue. Different, looser thresholds.
#
# Noise filter for discovery (per genotype SS=Shiraz / RR=Ruggeri):
#   L1 per-sample allele support: site depth >= DP_MIN, allele reads >= ALT_MIN,
#      VAF >= VAF_MIN, and (when alt reads >= 4) seen on BOTH strands.
#   L2 per-tissue reproducibility: allele supported (L1) in >= REPS_MIN replicates
#      of a tissue AND both strands present in the tissue's pooled reads. This is
#      what kills index-hopping / cross-sample contamination and one-off library
#      error (they do not reproduce across independent replicates).
#   L3 RNA-editing removal: drop an allele whose change is A->G or T->C and that
#      has NO genome (dipcall) support (A-to-I editing signature).
#   L4 diagnostic + recipient background: an allele in the donor's set but not the
#      recipient's is a marker only if, in the recipient self-graft, that tissue
#      has adequate coverage (pooled depth >= RECIP_DP_MIN in >= REPS_MIN reps)
#      and the allele is at background (fraction <= RECIP_BG).
#
# Evidence tier: genome+RNA (dipcall also carries the allele) vs RNA-only
# (assembly-gap / clone-specific — recovered only because of the RNA).
#
# Inputs:
#   <sites_dir>  : holds diagnostic_sites.tsv (Step 1 genome catalogue)
#   <cnt_dir>    : holds <sample>.cnt.tsv  (08; CHROM POS REF ALT AD ADF ADR)
#   <samplesheet>: metadata/step2_samples.tsv (sample group tissue rep ...)
#   [gff3]       : optional annotation for gene/exon tagging of RNA-only sites
# Params (positional, after gff): DP_MIN ALT_MIN VAF_MIN REPS_MIN RECIP_DP_MIN RECIP_BG
#   defaults:                       10     3       0.10    2        10           0.02
# Output: <sites_dir>/master_diagnostic_db_v3.tsv
# =====================================================================
suppressPackageStartupMessages(library(data.table))
setDTthreads(0)

a <- commandArgs(trailingOnly = TRUE)
sites_dir <- a[1]; cnt_dir <- a[2]; samplesheet <- a[3]
gffpath   <- if (length(a) >= 4 && nzchar(a[4]) && a[4] != "-") a[4] else NA_character_
DP_MIN       <- if (length(a) >= 5)  as.integer(a[5])  else 10L
ALT_MIN      <- if (length(a) >= 6)  as.integer(a[6])  else 3L
VAF_MIN      <- if (length(a) >= 7)  as.numeric(a[7])  else 0.10
REPS_MIN     <- if (length(a) >= 8)  as.integer(a[8])  else 2L
RECIP_DP_MIN <- if (length(a) >= 9)  as.integer(a[9])  else 10L
RECIP_BG     <- if (length(a) >= 10) as.numeric(a[10]) else 0.02

ss <- fread(samplesheet); setnames(ss, tolower(names(ss)))
ss <- ss[, .(sample, group, tissue)]

# ---- genome support: (CHROM,POS,allele) present as a genome allele ----
sites <- fread(file.path(sites_dir, "diagnostic_sites.tsv"))
gmelt <- function(col) {
  x <- sites[get(col) != "NA" & !is.na(get(col)), .(CHROM, POS, v = get(col))]
  if (!nrow(x)) return(x[, .(CHROM, POS, allele = character())])
  x[, .(allele = unlist(strsplit(v, ",", fixed = TRUE))), by = .(CHROM, POS)]
}
gsupp <- unique(rbind(gmelt("shiraz_alleles"), gmelt("ruggeri_alleles")))
gsupp[, gsup := TRUE]; setkey(gsupp, CHROM, POS, allele)

# ---- gene annotation: strand (for strand-aware editing) + id/name (for output) ----
genes_str <- NULL; genes_id <- NULL
if (!is.na(gffpath) && file.exists(gffpath)) {
  g <- fread(cmd = paste("grep -v '^#'", shQuote(gffpath)), header = FALSE, sep = "\t", quote = "")
  setnames(g, 1:9, c("chrom","src","type","start","end","score","strand","frame","attr"))
  g[, chrom := as.character(chrom)]
  gg <- g[type == "gene"]
  genes_str <- gg[, .(chrom, start, end, strand)]
  genes_id  <- gg[, .(chrom, start, end,
                      gene_id   = sub(".*ID=([^;]+).*", "\\1", attr),
                      gene_name = fifelse(grepl("Name=", attr),
                                          sub(".*Name=([^;]+).*", "\\1", attr), NA_character_))]
  setkey(genes_str, chrom, start, end); setkey(genes_id, chrom, start, end)
}

# ---- read all count tables, melt to per-allele rows (keep <*> for depth) ----
read_cnt <- function(f) {
  s <- sub("\\.cnt\\.tsv$", "", basename(f))
  d <- fread(f, header = FALSE, sep = "\t",
             col.names = c("CHROM","POS","REF","ALT","AD","ADF","ADR"),
             colClasses = list(character = c("CHROM","REF","ALT","AD","ADF","ADR")))
  if (!nrow(d)) return(NULL)
  d[, POS := as.integer(POS)]
  al <- strsplit(paste(d$REF, d$ALT, sep = ","), ",", fixed = TRUE)
  ad <- strsplit(d$AD,  ",", fixed = TRUE)
  af <- strsplit(d$ADF, ",", fixed = TRUE)
  ar <- strsplit(d$ADR, ",", fixed = TRUE)
  n  <- lengths(al)
  long <- data.table(CHROM = rep(d$CHROM, n), POS = rep(d$POS, n), REF = rep(d$REF, n),
                     allele = unlist(al),
                     ad  = as.integer(unlist(ad)),
                     adf = as.integer(unlist(af)),
                     adr = as.integer(unlist(ar)),
                     sample = s)
  long
}
cnt <- rbindlist(lapply(list.files(cnt_dir, "\\.cnt\\.tsv$", full.names = TRUE), read_cnt))
if (is.null(cnt) || !nrow(cnt)) stop("No *.cnt.tsv rows read from ", cnt_dir)
cnt <- merge(cnt, ss, by = "sample")

# site depth per (sample, position) = sum over all alleles (incl <*>)
cnt[, dp := sum(ad), by = .(sample, CHROM, POS)]

# ---- L1: per-sample allele support (real alleles only) ----
cnt[, is_star := allele == "<*>"]
cnt[, L1 := !is_star & dp >= DP_MIN & ad >= ALT_MIN & (ad / dp) >= VAF_MIN &
            (ad < 4 | (adf >= 1 & adr >= 1))]

# ---- L2: per (group,tissue,pos,allele) reproducibility + pooled strand ----
L2 <- cnt[is_star == FALSE,
          .(nsupp = sum(L1), padf = sum(adf), padr = sum(adr),
            pooled_cnt = sum(ad), REF = REF[1]),
          by = .(group, tissue, CHROM, POS, allele)]
L2[, present := nsupp >= REPS_MIN & padf >= 1 & padr >= 1]

# ---- strand per position (editing direction); conservative (both) if no gene ----
L2[, `:=`(has_plus = TRUE, has_minus = TRUE)]
if (!is.null(genes_str)) {
  upos <- unique(L2[, .(CHROM, POS)])
  dq <- data.table(chrom = upos$CHROM, start = upos$POS, end = upos$POS)
  setkey(dq, chrom, start, end)
  ov <- foverlaps(dq, genes_str, type = "any", nomatch = NULL)
  stab <- ov[, .(hp = any(strand == "+"), hm = any(strand == "-")),
             by = .(CHROM = chrom, POS = i.start)]
  L2[stab, on = .(CHROM, POS), `:=`(has_plus = i.hp, has_minus = i.hm)]
}

# ---- L3: remove RNA-editing alleles (strand-aware A->G / T->C) without genome support ----
# A-to-I editing shows as A->G on (+)-strand genes and T->C on (-)-strand genes
# in reference coordinates. This uses the GENE's strand (from the gff), which is
# independent of whether the RNA-seq library is strand-specific.
L2 <- merge(L2, gsupp, by = c("CHROM","POS","allele"), all.x = TRUE)
L2[is.na(gsup), gsup := FALSE]
L2[present == TRUE & gsup == FALSE &
     ((has_plus  & REF == "A" & allele == "G") |
      (has_minus & REF == "T" & allele == "C")),
   present := FALSE]

# ---- pooled depth & replicate coverage per (group,tissue,pos) ----
dptab <- unique(cnt[, .(sample, group, tissue, CHROM, POS, dp)])
cov <- dptab[, .(pool_dp = sum(dp), nrep_cov = sum(dp >= DP_MIN)),
             by = .(group, tissue, CHROM, POS)]
setkey(cov, group, tissue, CHROM, POS)

# allele sets that passed L1-L3 (carry REF through for the final output)
pres <- L2[present == TRUE, .(group, tissue, CHROM, POS, REF, allele, pooled_cnt, gsup)]

# ---- L4 (v3): donor allele from ANY donor tissue, checked per recipient tissue ----
dirs <- data.table(direction = c("shiraz","ruggeri"),
                    donor     = c("SS","RR"),
                    recip     = c("RR","SS"))
res <- list()
for (i in 1:2) {
  dn <- dirs$donor[i]; rc <- dirs$recip[i]; dr <- dirs$direction[i]
  # donor allele present (L1-L3) in at least one donor tissue
  dps <- pres[group == dn,
              .(REF = REF[1], pooled_cnt = sum(pooled_cnt), gsup = any(gsup),
                donor_tissues = paste(sort(unique(tissue)), collapse = ",")),
              by = .(CHROM, POS, allele)]
  # test it in every tissue of the recipient
  rtis <- sort(unique(ss[group == rc, tissue]))
  cand <- dps[rep(seq_len(nrow(dps)), each = length(rtis))]
  cand[, tissue := rep(rtis, times = nrow(dps))]
  rps <- pres[group == rc, .(tissue, CHROM, POS, allele)]    # recipient present alleles
  rps[, rpresent := TRUE]
  cand <- merge(cand, rps, by = c("tissue","CHROM","POS","allele"), all.x = TRUE)
  cand <- cand[is.na(rpresent)]                              # allele NOT in recipient set
  # recipient coverage + background (same tissue)
  rcov <- cov[group == rc, .(tissue, CHROM, POS, recip_dp = pool_dp, recip_nrep = nrep_cov)]
  cand <- merge(cand, rcov, by = c("tissue","CHROM","POS"), all.x = TRUE)
  rcnt <- cnt[group == rc & is_star == FALSE,
              .(recip_cnt = sum(ad)), by = .(tissue, CHROM, POS, allele)]
  cand <- merge(cand, rcnt, by = c("tissue","CHROM","POS","allele"), all.x = TRUE)
  cand[is.na(recip_dp), recip_dp := 0L][is.na(recip_nrep), recip_nrep := 0L]
  cand[is.na(recip_cnt), recip_cnt := 0L]
  cand[, recip_frac := fifelse(recip_dp > 0, recip_cnt / recip_dp, NA_real_)]
  cand <- cand[recip_dp >= RECIP_DP_MIN & recip_nrep >= REPS_MIN &
               !is.na(recip_frac) & recip_frac <= RECIP_BG]
  cand[, `:=`(direction = dr, donor_cnt = pooled_cnt,
              tier = fifelse(gsup, "genome+RNA", "RNA-only"))]
  res[[i]] <- cand[, .(CHROM, POS, REF, tissue, direction, allele,
                       donor_cnt, recip_dp, recip_frac = round(recip_frac, 4), tier,
                       donor_tissues)]
}
db <- rbindlist(res)

# ---- gene id / name tag (from the genes loaded above) ----
if (!is.null(genes_id) && nrow(db)) {
  # FIX (v3b): setkey() re-sorts dq, so foverlaps() returns rows in POSITION order,
  # not in db row order. Carry a row id and restore db order before assigning.
  # (Earlier versions assigned gene ids to the wrong rows for RNA-only markers.)
  dq <- data.table(chrom = db$CHROM, start = db$POS, end = db$POS, rid = seq_len(nrow(db)))
  setkey(dq, chrom, start, end)
  gi <- foverlaps(dq, genes_id, type = "any", mult = "first")
  setorder(gi, rid)
  stopifnot(nrow(gi) == nrow(db), identical(gi$rid, seq_len(nrow(db))))
  db[, gene_id := gi$gene_id][, gene_name := gi$gene_name]
} else db[, `:=`(gene_id = NA_character_, gene_name = NA_character_)]

setorder(db, CHROM, POS, tissue, direction)
fwrite(db, file.path(sites_dir, "master_diagnostic_db_v3.tsv"), sep = "\t", na = "NA", quote = FALSE)

message("Master diagnostic markers (site x tissue x direction x allele): ", nrow(db))
message("  shiraz-private: ", db[direction=="shiraz", .N],
        " ; ruggeri-private: ", db[direction=="ruggeri", .N])
message("  tier genome+RNA: ", db[tier=="genome+RNA", .N],
        " ; RNA-only (gap/clone): ", db[tier=="RNA-only", .N])
message("  distinct positions: ", uniqueN(db[, .(CHROM, POS)]),
        " ; tissues: ", paste(sort(unique(db$tissue)), collapse=","))
message("  -> ", file.path(sites_dir, "master_diagnostic_db_v3.tsv"))
