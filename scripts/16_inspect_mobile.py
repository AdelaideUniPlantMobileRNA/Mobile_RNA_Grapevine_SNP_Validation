#!/usr/bin/env python3
# =====================================================================
# Step 4c: quick checks of the mobile calls (no new Bayes factors computed).
#  (1) GLOBAL check (Paajanen et al. 2025 recommendation): in each heterograft
#      sample, fraction of reads carrying the mobile allele over all scored markers,
#      compared with the same fraction in the local self-graft at the same markers.
#      If heterograft ~= self-graft, the sample as a whole shows no excess.
#  (2) PER-GENE check for every called mobile gene, per replicate:
#      SNPs covered, SNPs with mobile reads, mobile-allele fraction in heterograft
#      vs local self-graft, and how much of the gene's evidence comes from ONE SNP
#      (top_snp_share close to 1 = driven by a single SNP -> check for mismapping).
#  (3) overlap with an older call list (optional).
# Usage: python 16_inspect_mobile.py <out_dir> [old_calls.tsv|-] [exclude=SRR4] [hom_cap=300]
#        (use hom_cap 2000 for results made by 13 v3)
# =====================================================================
import sys, os, glob
import numpy as np
import pandas as pd

OUT = sys.argv[1]
OLD = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] != "-" else None
EXCLUDE = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else {"SRR4"}
CAP = int(sys.argv[4]) if len(sys.argv) > 4 else 300   # = self-graft cap used in 13

calls = pd.read_csv(os.path.join(OUT, "mobile_transcripts.tsv"), sep="\t")
mob = calls[calls["mobile"]][["graft", "tissue", "direction", "gene_id"]]
keys = set(map(tuple, mob.values))

glob_rows, gene_rows = [], []
for f in sorted(glob.glob(os.path.join(OUT, "per_sample", "*.snp.tsv"))):
    s = os.path.basename(f).replace(".snp.tsv", "")
    if s in EXCLUDE:
        continue
    d = pd.read_csv(f, sep="\t", dtype={"CHROM": str})
    # (1) global fraction, heterograft depth capped like the self-graft
    Nc = np.minimum(d["N"], CAP); nc = d["n"] * Nc / d["N"].clip(lower=1)
    hom_ok = d["Nh1"] > 0
    glob_rows.append(dict(
        sample=s, graft=d["graft"].iat[0], tissue=d["tissue"].iat[0],
        direction=d["direction"].iat[0], rep=d["rep"].iat[0], n_snp=len(d),
        het_mobile_frac=nc.sum() / Nc.sum(),
        hom_mobile_frac=d.loc[hom_ok, "nh1"].sum() / d.loc[hom_ok, "Nh1"].sum(),
        het_snps_with_mobile_reads=(d["n"] > 0).mean(),
        hom_snps_with_mobile_reads=(d.loc[hom_ok, "nh1"] > 0).mean()))
    # (2) per-gene detail for called genes
    k = d[["graft", "tissue", "direction", "gene_id"]].apply(tuple, axis=1).isin(keys)
    for gid, x in d[k].groupby("gene_id"):
        pos = x["log10BF"].clip(lower=0)
        gene_rows.append(dict(
            graft=x["graft"].iat[0], tissue=x["tissue"].iat[0],
            direction=x["direction"].iat[0], gene_id=gid, sample=s, rep=x["rep"].iat[0],
            n_snp=len(x), snp_with_mobile_reads=int((x["n"] > 0).sum()),
            snp_BF_ge1=int((x["log10BF"] >= 1).sum()),
            het_frac=round(x["n"].sum() / x["N"].sum(), 4),
            hom_frac=round(x["nh1"].sum() / max(x["Nh1"].sum(), 1), 4),
            het_N=int(x["N"].sum()), het_n=int(x["n"].sum()),
            sum_log10BF=round(x["log10BF"].sum(), 2),
            top_snp_share=round(pos.max() / pos.sum(), 2) if pos.sum() > 0 else np.nan))

g = pd.DataFrame(glob_rows).sort_values(["graft", "tissue", "rep"])
g["ratio_het_vs_hom"] = (g["het_mobile_frac"] / g["hom_mobile_frac"]).round(2)
g.to_csv(os.path.join(OUT, "check_global_fraction.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print("==== (1) GLOBAL mobile-allele fraction: heterograft vs local self-graft ====")
print(g[["sample", "graft", "tissue", "direction", "n_snp", "het_mobile_frac",
         "hom_mobile_frac", "ratio_het_vs_hom"]].to_string(index=False,
         float_format=lambda v: f"{v:.5f}"))

gd = pd.DataFrame(gene_rows).sort_values(["graft", "tissue", "direction", "gene_id", "rep"])
gd.to_csv(os.path.join(OUT, "check_mobile_genes.tsv"), sep="\t", index=False)
print("\n==== (2) called mobile genes, per replicate ====")
print(gd.to_string(index=False))

summ = gd[gd["sum_log10BF"] >= 1].groupby(["graft", "tissue", "direction", "gene_id"]).agg(
    reps=("rep", "nunique"), med_snps=("n_snp", "median"),
    med_snps_mobile=("snp_with_mobile_reads", "median"),
    med_het_frac=("het_frac", "median"), med_hom_frac=("hom_frac", "median"),
    max_top_snp_share=("top_snp_share", "max")).reset_index()
print("\n==== summary over passing replicates ====")
print(summ.to_string(index=False))
print("\n  single-SNP driven (top_snp_share >= 0.8 in a passing rep):",
      int((summ["max_top_snp_share"] >= 0.8).sum()), "of", len(summ))
print("  by chromosome:", summ["gene_id"].str.extract(r"(chr\d+)")[0].value_counts().to_dict())

if OLD and os.path.exists(OLD):
    o = pd.read_csv(OLD, sep="\t"); o = o[o["mobile"]]
    ok = set(map(tuple, o[["graft", "tissue", "direction", "gene_id"]].values))
    print(f"\n==== (3) overlap with old calls: {len(ok & keys)} of {len(ok)} old calls are in the new {len(keys)}")
    for t in sorted(ok): print("   ", t, "IN NEW" if t in keys else "not in new")
print("\n  -> check_global_fraction.tsv, check_mobile_genes.tsv in", OUT)
