#!/usr/bin/env python3
# =====================================================================
# Step 4a (v3): baymobil Bayes factors for ONE heterograft sample.
#
# v3 = v2 plus:
#   * RNA-only markers come from master_diagnostic_db_v3.tsv (09 v3: donor allele
#     confirmed in ANY donor tissue, e.g. RR root, not only the same tissue).
#   * "SSL shows it -> remove" now applies at ANY self-graft depth
#     (local self-graft has >= 1 read of the allele AND fraction > 2%).
#     v2 required depth >= 10, so low-depth SSL sites with e.g. 8/17 Ruggeri reads
#     slipped through and gave false calls.
#   * self-graft depth cap raised 300 -> 2000 (error rate estimated more precisely,
#     so low mobile fractions in highly expressed genes can be detected).
#   * TWO results per sample, reported side by side later (13b v3):
#       <sample>.snp.tsv  : per-SNP Bayes factors (baymobil original; summed per gene in 13b)
#       <sample>.gene.tsv : gene-pooled Bayes factor - counts summed over SNPs that are
#                           >= THIN bp apart (so one read pair is not counted twice),
#                           then ONE baymobil Bayes factor per gene.
#
# Marker set for this sample (tissue T, mobile genotype D):
#   (a) genome markers : Step 1 diagnostic_sites.tsv, exonic SNPs, allele private to D
#   (b) RNA-only markers: master_diagnostic_db_v3 rows for (T, D), tier RNA-only
#   tested if the heterograft sample has >= 1 read there; removed if the LOCAL
#   self-graft of tissue T shows the allele (see above).
#
# baymobil inputs per SNP (baymobil's own fasterpostN2, maths unchanged):
#   N , n   = heterograft total / mobile-allele reads (this sample)
#   Nh1,nh1 = LOCAL self-graft, same tissue (e.g. SSL reps): total / mobile-allele reads
#   Nh2,nh2 = DISTAL self-graft, all tissues pooled (e.g. all RR): total / non-mobile reads
#   no self-graft reads -> Nh=nh=0 -> flat Beta(1,1) prior (baymobil default behaviour)
#
# Usage:
#   python 13_detect_mobile_v3.py <master_db_v3> <diag_sites> <hom_cnt_dir> <het_cnt_dir> \
#          <hom_sheet> <het_sheet> <out_dir> <sample> [nmax=10] [hom_cap=2000] [thin_bp=300]
# =====================================================================
import sys, os
import numpy as np
import pandas as pd
from multiprocessing import Pool

try:                                   # baymobil's core Bayes-factor function
    from baymobil.baymobil import fasterpostN2
except Exception:
    import baymobil
    fasterpostN2 = getattr(baymobil, "fasterpostN2")

(master_db, diag_sites, HOMDIR, HETDIR, hom_sheet, het_sheet, OUT, SAMPLE) = sys.argv[1:9]
NMAX     = int(sys.argv[9])  if len(sys.argv) > 9  else 10
HOM_CAP  = int(sys.argv[10]) if len(sys.argv) > 10 else 2000
THIN     = int(sys.argv[11]) if len(sys.argv) > 11 else 300
RM_FRAC  = 0.02                      # remove marker if local self-graft shows the allele
POOL_N_CAP = 5000                    # gene-pooled heterograft depth cap (speed only)
THREADS  = int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))

outdir = os.path.join(OUT, "per_sample"); os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, f"{SAMPLE}.snp.tsv")
genefile = os.path.join(outdir, f"{SAMPLE}.gene.tsv")
if os.path.exists(outfile) and os.path.exists(genefile):
    print(f">> {SAMPLE}: already done, skipping"); sys.exit(0)

hom = pd.read_csv(hom_sheet, sep="\t")
het = pd.read_csv(het_sheet, sep="\t")
hs = het[het["sample"] == SAMPLE]
if hs.empty:
    sys.exit(f"sample {SAMPLE} not in {het_sheet}")
hs = hs.iloc[0]
T, direction = hs["tissue"], hs["direction"]
local_geno, distal_geno = hs["local_geno"], hs["mobile_geno"]

def is_true(s):
    return s.astype(str).str.upper().isin(["TRUE", "T", "1"])

# ---- (a) genome markers: exonic SNP alleles private to the mobile genotype ----
pcol = "ruggeri_private" if direction == "ruggeri" else "shiraz_private"
g = pd.read_csv(diag_sites, sep="\t", dtype={"CHROM": str},
                usecols=["CHROM", "POS", "REF", pcol, "gene_id", "is_exonic"])
g = g[is_true(g["is_exonic"]) & g[pcol].notna() & g["gene_id"].notna()]
g["allele"] = g[pcol].astype(str).str.split(",")
g = g.explode("allele")
g = g[(g["REF"].astype(str).str.len() == 1) & (g["allele"].str.len() == 1)]
g = g[["CHROM", "POS", "allele", "gene_id"]].assign(source="genome")

# ---- (b) RNA-only markers from the Step 2 master DB (this tissue + direction) ----
m = pd.read_csv(master_db, sep="\t", dtype={"CHROM": str},
                usecols=["CHROM", "POS", "REF", "tissue", "direction", "allele", "tier", "gene_id"])
m = m[(m["tissue"] == T) & (m["direction"] == direction) & (m["tier"] == "RNA-only")]
m = m[(m["REF"].astype(str).str.len() == 1) & (m["allele"].astype(str).str.len() == 1)]
m = m[["CHROM", "POS", "allele", "gene_id"]].assign(source="RNA-only")

mk = pd.concat([g, m], ignore_index=True)
mk["POS"] = mk["POS"].astype(int)
mk = mk.drop_duplicates(["CHROM", "POS", "allele"])
keep_pos = mk[["CHROM", "POS"]].drop_duplicates()
print(f">> {SAMPLE} ({T}/{direction}, local={local_geno}, distal={distal_geno}): "
      f"{len(mk)} SNP markers (genome {int((mk.source=='genome').sum())}, "
      f"RNA-only {int((mk.source=='RNA-only').sum())})", flush=True)

def read_cnt_long(path):
    d = pd.read_csv(path, sep="\t", header=None,
                    names=["CHROM", "POS", "REF", "ALT", "AD", "ADF", "ADR"],
                    dtype={0: str}, usecols=[0, 1, 2, 3, 4])
    d["POS"] = d["POS"].astype(int)
    d = d.merge(keep_pos, on=["CHROM", "POS"], how="inner")
    # SNP records only (skip mpileup's extra indel record at the same position)
    alts = d["ALT"].astype(str).str.split(",")
    snp_line = (d["REF"].astype(str).str.len() == 1) & \
               alts.apply(lambda a: all(len(x) == 1 or x == "<*>" for x in a))
    d = d[snp_line]
    d["allele"] = (d["REF"].astype(str) + "," + d["ALT"].astype(str)).str.split(",")
    d["cnt"]    = d["AD"].astype(str).str.split(",")
    long = d[["CHROM", "POS", "allele", "cnt"]].explode(["allele", "cnt"])
    long["cnt"] = long["cnt"].astype(int)
    tot  = long.groupby(["CHROM", "POS"], as_index=False)["cnt"].sum().rename(columns={"cnt": "tot"})
    long = long.groupby(["CHROM", "POS", "allele"], as_index=False)["cnt"].sum()
    return long, tot

def pooled_hom(geno, tissue=None):
    sel = hom["group"] == geno
    if tissue is not None:
        sel &= hom["tissue"] == tissue
    L, Tt = [], []
    for s in hom.loc[sel, "sample"]:
        p = os.path.join(HOMDIR, f"{s}.cnt.tsv")
        if os.path.exists(p):
            a, b = read_cnt_long(p); L.append(a); Tt.append(b)
        else:
            print(f"   warning: missing {p}", flush=True)
    L  = pd.concat(L).groupby(["CHROM", "POS", "allele"], as_index=False)["cnt"].sum()
    Tt = pd.concat(Tt).groupby(["CHROM", "POS"], as_index=False)["tot"].sum()
    return L, Tt

hl, ht     = read_cnt_long(os.path.join(HETDIR, f"{SAMPLE}.cnt.tsv"))
locL, locT = pooled_hom(local_geno, T)        # e.g. SSL reps (same tissue)
disL, disT = pooled_hom(distal_geno)          # e.g. all RR samples, all tissues

df = mk.merge(hl.rename(columns={"cnt": "n"}),      on=["CHROM", "POS", "allele"], how="left") \
       .merge(ht.rename(columns={"tot": "N"}),      on=["CHROM", "POS"],           how="left") \
       .merge(locL.rename(columns={"cnt": "nh1"}),  on=["CHROM", "POS", "allele"], how="left") \
       .merge(locT.rename(columns={"tot": "Nh1"}),  on=["CHROM", "POS"],           how="left") \
       .merge(disL.rename(columns={"cnt": "dcnt"}), on=["CHROM", "POS", "allele"], how="left") \
       .merge(disT.rename(columns={"tot": "Nh2"}),  on=["CHROM", "POS"],           how="left")
for c in ["n", "N", "nh1", "Nh1", "dcnt", "Nh2"]:
    df[c] = df[c].fillna(0).astype(int)
df["nh2"] = df["Nh2"] - df["dcnt"]

# only markers expressed in THIS heterograft sample
n_all = len(df)
df = df[df["N"] >= 1].copy()
n_expr = len(df)
# remove markers where the local self-graft (same tissue) shows the allele
shows = (df["nh1"] >= 1) & (df["nh1"] > RM_FRAC * df["Nh1"])
df = df[~shows].copy()
print(f"   markers: {n_all} total -> {n_expr} with >=1 read in {SAMPLE} -> "
      f"{len(df)} after removing {int(shows.sum())} seen in local self-graft (>= 1 read and > {RM_FRAC:.0%})",
      flush=True)
print(f"   of these: no local self-graft reads (Nh1=0): {int((df.Nh1 == 0).sum())}; "
      f"no distal self-graft reads (Nh2=0): {int((df.Nh2 == 0).sum())}", flush=True)

# ---- gene-pooled counts (uncapped here; capped just before the Bayes factor) ----
srt = df.sort_values(["gene_id", "CHROM", "POS"]).drop_duplicates(["gene_id", "CHROM", "POS"])
pos_arr = srt["POS"].to_numpy(); keep = np.zeros(len(srt), dtype=bool)
for idx in srt.groupby(["gene_id", "CHROM"], sort=False).indices.values():
    last = -10**12                          # greedy by position; independent of n
    for j in idx:
        if pos_arr[j] - last >= THIN:
            keep[j] = True; last = pos_arr[j]
thin = srt[keep]
gp = thin.groupby("gene_id", as_index=False).agg(
    n_snp_used=("POS", "size"), N=("N", "sum"), n=("n", "sum"),
    Nh1=("Nh1", "sum"), nh1=("nh1", "sum"), Nh2=("Nh2", "sum"), nh2=("nh2", "sum"))
gp["n_snp_all"] = gp["gene_id"].map(df.groupby("gene_id").size())

# cap self-graft depth (cost of logbeta grows with it; error rate is well estimated
# from ~300 reads). Heterograft N is NOT capped.
for Ncol, ncol in (("Nh1", "nh1"), ("Nh2", "nh2")):
    over = df[Ncol] > HOM_CAP
    df.loc[over, ncol] = np.round(df.loc[over, ncol] * HOM_CAP / df.loc[over, Ncol]).astype(int)
    df.loc[over, Ncol] = HOM_CAP

cols = ["Nh1", "nh1", "Nh2", "nh2", "N", "n"]
uniq = df[cols].drop_duplicates().astype(int)
tuples = [tuple(int(v) for v in r) + (NMAX,) for r in uniq.itertuples(index=False)]
print(f"   {len(df)} SNPs scored, {len(tuples)} unique input tuples, {THREADS} workers", flush=True)

def _bf(t):
    meanN2, N2max, lbf = fasterpostN2(*t)
    return (float(meanN2), int(N2max), float(lbf))

with Pool(THREADS) as pool:
    res = pool.map(_bf, tuples, chunksize=max(1, len(tuples) // (THREADS * 20)))

uniq[["meanN2", "N2max", "log10BF"]] = pd.DataFrame(res, index=uniq.index)
df = df.merge(uniq, on=cols, how="left")
df["sample"] = SAMPLE; df["graft"] = hs["graft"]; df["tissue"] = T
df["direction"] = direction; df["rep"] = hs["rep"]
df = df[["sample", "graft", "tissue", "direction", "rep", "CHROM", "POS", "allele", "gene_id",
         "source", "N", "n", "Nh1", "nh1", "Nh2", "nh2", "meanN2", "N2max", "log10BF"]]
tmp = outfile + ".tmp"
df.to_csv(tmp, sep="\t", index=False); os.replace(tmp, outfile)   # atomic: no half files

# ---- gene-pooled Bayes factor (one per gene) ----
def capcol(d, Ncol, ncol, cap):
    over = d[Ncol] > cap
    d.loc[over, ncol] = np.round(d.loc[over, ncol] * cap / d.loc[over, Ncol]).astype(int)
    d.loc[over, Ncol] = cap
capcol(gp, "Nh1", "nh1", HOM_CAP); capcol(gp, "Nh2", "nh2", HOM_CAP)
gp["N_raw"], gp["n_raw"] = gp["N"], gp["n"]
capcol(gp, "N", "n", POOL_N_CAP)
gtup = [tuple(int(v) for v in r) + (NMAX,) for r in gp[cols].itertuples(index=False)]
print(f"   gene-pooled: {len(gp)} genes, SNPs >= {THIN} bp apart", flush=True)
with Pool(THREADS) as pool:
    gres = pool.map(_bf, gtup, chunksize=max(1, len(gtup) // (THREADS * 20)))
gp[["meanN2", "N2max", "log10BF"]] = pd.DataFrame(gres, index=gp.index)
gp["sample"] = SAMPLE; gp["graft"] = hs["graft"]; gp["tissue"] = T
gp["direction"] = direction; gp["rep"] = hs["rep"]
gp = gp[["sample", "graft", "tissue", "direction", "rep", "gene_id", "n_snp_all", "n_snp_used",
         "N_raw", "n_raw", "N", "n", "Nh1", "nh1", "Nh2", "nh2", "meanN2", "N2max", "log10BF"]]
tmp = genefile + ".tmp"
gp.to_csv(tmp, sep="\t", index=False); os.replace(tmp, genefile)
print(f">> {SAMPLE}: gene-pooled done, {(gp['log10BF'] >= 1).sum()} genes with log10BF>=1 -> {genefile}",
      flush=True)
print(f">> {SAMPLE}: done, {df['gene_id'].nunique()} genes tested, "
      f"{(df['log10BF'] >= 1).sum()} SNPs with log10BF>=1 -> {outfile}", flush=True)
