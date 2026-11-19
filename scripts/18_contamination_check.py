#!/usr/bin/env python3
# =====================================================================
# Step 4d: contamination / cross-sample check (after Paajanen et al. 2025).
#
# Idea: if donor-allele reads in a heterograft sample came from CONTAMINATION
# (tissue carry-over from the graft partner, or index hopping from donor
# libraries), their number at each marker would be PROPORTIONAL to how much the
# donor source expresses that marker - across ALL genes. Real mobility affects a
# few genes only.
#
# For every heterograft sample, at its tested SNP markers, we fit
#     n_i  ~  Poisson( a * N_i  +  c * D_i )
#   n_i = donor-allele reads in the heterograft sample
#   N_i = total reads in the heterograft sample   (a = sequencing-error rate)
#   D_i = donor-allele reads in a candidate contamination SOURCE (c = contamination rate)
# Sources:
#   partner : the other side of the SAME graft type (e.g. for SRL: SR root samples), pooled
#   donorHG : donor-genotype homograft of the SAME tissue (e.g. for SRL: RR leaf), pooled
# Fitted on all markers EXCEPT the called mobile genes. Output:
#   contam_share = c*sum(D) / (a*sum(N) + c*sum(D))  = share of donor-allele reads
#                  explained by that source.
# NULL: the same model fitted on each recipient HOMOGRAFT replicate of the same tissue
#   (e.g. SSL1..SSL4) at the same markers. A homograft has no graft partner, so its
#   contam_share is the background level (index hopping / shared batch effects).
#   Heterograft contam_share >> homograft  -> contamination from that source.
# Called genes: observed donor reads vs reads expected from error + contamination.
#
# Usage: python 18_contamination_check.py <out_dir> <hom_cnt_dir> <het_cnt_dir> \
#          <hom_sheet> <het_sheet> [exclude=SRR4]
# =====================================================================
import sys, os
import numpy as np
import pandas as pd

OUT, HOMDIR, HETDIR, hom_sheet, het_sheet = sys.argv[1:6]
EXCLUDE = set(sys.argv[6].split(",")) if len(sys.argv) > 6 else {"SRR4"}
hom = pd.read_csv(hom_sheet, sep="\t")
het = pd.read_csv(het_sheet, sep="\t")
het = het[~het["sample"].isin(EXCLUDE)]
calls = pd.read_csv(os.path.join(OUT, "mobile_transcripts.tsv"), sep="\t")
calls = calls[calls["mobile"]]
BASES = ["A", "C", "G", "T"]

_cache = {}
def read_wide(path):
    """per position: reads of A,C,G,T and total depth (SNP records only)."""
    if path in _cache:
        return _cache[path]
    d = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 3, 4],
                    names=["CHROM", "POS", "REF", "ALT", "AD"], dtype={0: str})
    alts = d["ALT"].astype(str).str.split(",")
    keep = (d["REF"].astype(str).str.len() == 1) & \
           alts.apply(lambda a: all(len(x) == 1 or x == "<*>" for x in a))
    d = d[keep]
    al = (d["REF"] + "," + d["ALT"]).str.split(",")
    ad = d["AD"].astype(str).str.split(",")
    long = pd.DataFrame({"CHROM": np.repeat(d["CHROM"].values, al.str.len()),
                         "POS": np.repeat(d["POS"].values, al.str.len()),
                         "allele": np.concatenate(al.values),
                         "cnt": np.concatenate(ad.values).astype(np.int32)})
    w = long.pivot_table(index=["CHROM", "POS"], columns="allele", values="cnt",
                         aggfunc="sum", fill_value=0)
    for b in BASES:
        if b not in w.columns:
            w[b] = 0
    w["DP"] = w.sum(axis=1)
    w = w[BASES + ["DP"]].astype(np.int32)
    _cache[path] = w
    return w

def counts_at(paths, mk):
    """pooled donor-allele reads and depth over `paths` at markers mk."""
    idx = pd.MultiIndex.from_frame(mk[["CHROM", "POS"]])
    don = np.zeros(len(mk), dtype=np.int64); dp = np.zeros(len(mk), dtype=np.int64)
    for p in paths:
        if not os.path.exists(p):
            continue
        w = read_wide(p).reindex(idx, fill_value=0)
        dp += w["DP"].values
        don += np.choose(mk["bidx"].values, [w[b].values for b in BASES])
    return don, dp

def fit(n, N, D, iters=300):
    """ML-EM for n ~ Poisson(a*N + c*D), a,c >= 0."""
    a = max(n.sum(), 1) / max(N.sum(), 1) / 2; c = a
    sN, sD = N.sum(), max(D.sum(), 1)
    for _ in range(iters):
        mu = a * N + c * D
        ok = mu > 0
        r = np.zeros_like(mu); r[ok] = n[ok] / mu[ok]
        a = a * (r * N).sum() / sN
        c = c * (r * D).sum() / sD
    share = c * D.sum() / max(a * N.sum() + c * D.sum(), 1e-12)
    return a, c, share

res, gene_rows = [], []
for _, s in het.iterrows():
    f = os.path.join(OUT, "per_sample", f"{s['sample']}.snp.tsv")
    mk = pd.read_csv(f, sep="\t", dtype={"CHROM": str},
                     usecols=["CHROM", "POS", "allele", "gene_id", "N", "n"])
    mk = mk[mk["allele"].isin(BASES)].reset_index(drop=True)
    mk["bidx"] = mk["allele"].map({b: i for i, b in enumerate(BASES)})
    called = set(calls[(calls.graft == s.graft) & (calls.tissue == s.tissue) &
                       (calls.direction == s.direction)]["gene_id"])
    bg = ~mk["gene_id"].isin(called)
    src = {
        "partner": [os.path.join(HETDIR, f"{x}.cnt.tsv") for x in
                    het[(het.graft == s.graft) & (het.local_geno == s.mobile_geno)]["sample"]],
        "donorHG": [os.path.join(HOMDIR, f"{x}.cnt.tsv") for x in
                    hom[(hom.group == s.mobile_geno) & (hom.tissue == s.tissue)]["sample"]],
    }
    null_reps = hom[(hom.group == s.local_geno) & (hom.tissue == s.tissue)]["sample"].tolist()
    n = mk["n"].values.astype(float); N = mk["N"].values.astype(float)
    for sname, paths in src.items():
        D, _ = counts_at(paths, mk); D = D.astype(float)
        a, c, share = fit(n[bg], N[bg], D[bg])
        nulls = []
        for r in null_reps:
            hn, hN = counts_at([os.path.join(HOMDIR, f"{r}.cnt.tsv")], mk)
            m = hN > 0
            nulls.append(fit(hn[m & bg].astype(float), hN[m & bg].astype(float), D[m & bg])[2])
        res.append(dict(sample=s["sample"], graft=s.graft, tissue=s.tissue,
                        direction=s.direction, source=sname, n_src_files=len(paths),
                        error_rate=a, contam_rate=c, contam_share=share,
                        homograft_share_median=np.median(nulls) if nulls else np.nan,
                        homograft_share_max=np.max(nulls) if nulls else np.nan))
        # called genes in this sample: observed vs expected donor reads
        for gid in called:
            g = mk["gene_id"] == gid
            if not g.any():
                continue
            exp = a * N[g].sum() + c * D[g].sum()
            gene_rows.append(dict(sample=s["sample"], graft=s.graft, tissue=s.tissue,
                                  direction=s.direction, gene_id=gid, source=sname,
                                  obs_n=int(n[g].sum()), expected_n=round(exp, 2),
                                  obs_over_exp=round(n[g].sum() / max(exp, 1e-9), 1)))
    print(f">> {s['sample']} done", flush=True)

R = pd.DataFrame(res)
R.to_csv(os.path.join(OUT, "check_contamination.tsv"), sep="\t", index=False)
G = pd.DataFrame(gene_rows)
G.to_csv(os.path.join(OUT, "check_contamination_genes.tsv"), sep="\t", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print("\n==== share of donor-allele reads explained by each contamination source ====")
print("     (heterograft vs the same model on recipient homograft replicates = background)")
print(R[["sample", "graft", "tissue", "direction", "source", "error_rate", "contam_share",
         "homograft_share_median", "homograft_share_max"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4f}"))
if len(G):
    Gp = G.groupby(["graft", "tissue", "direction", "gene_id", "source"])["obs_over_exp"].median().unstack()
    print("\n==== called mobile genes: observed / expected(error + contamination), median over samples ====")
    for col in Gp.columns:
        v = Gp[col].dropna()
        print(f"  source {col}: {len(v)} genes; obs/exp < 2: {(v < 2).sum()} | 2-5: "
              f"{((v >= 2) & (v < 5)).sum()} | >= 5: {(v >= 5).sum()}")
print("\n  -> check_contamination.tsv, check_contamination_genes.tsv in", OUT)
