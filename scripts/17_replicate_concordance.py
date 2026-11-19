#!/usr/bin/env python3
# =====================================================================
# Replicate-concordance check (no new Bayes factors computed).
# Question: are the "mobile" genes (pass in >= 2 replicates) MORE reproducible
# than expected if each replicate passed genes at random?
#   For each sample s: p_s = fraction of its tested genes with gene log10BF >= 1.
#   If passes were independent noise, a gene tested in reps {s1..sk} passes in
#   >= m reps with probability P(>= m) under a Poisson-binomial(p_s1..p_sk).
#   Summing over genes gives the EXPECTED number of calls by chance.
# Observed >> expected  -> replicate agreement is real (mobility, or a consistent
#                          artefact such as paralog mis-mapping - not proof of mobility).
# Observed ~= expected  -> calls are indistinguishable from independent noise.
# Usage: python 17_replicate_concordance.py <out_dir>
# =====================================================================
import sys, os
import numpy as np
import pandas as pd

OUT = sys.argv[1]
KEY = ["graft", "tissue", "direction"]

def poibin_tail(ps):
    """P(X >= m) for m = 0..k, X = sum of independent Bernoulli(ps)."""
    dist = np.array([1.0])
    for p in ps:
        dist = np.convolve(dist, [1 - p, p])
    return dist[::-1].cumsum()[::-1]          # tail[m] = P(X >= m)

def run(per_rep, score_col, label):
    d = per_rep.copy()
    d["pass"] = d[score_col] >= 1
    p_s = d.groupby("sample")["pass"].mean()
    d["p"] = d["sample"].map(p_s)
    rows = []
    for combo, g in d.groupby(KEY):
        exp2 = exp3 = exp4 = 0.0
        obs = g.groupby("gene_id")["pass"].sum()
        for gid, x in g.groupby("gene_id"):
            t = poibin_tail(x["p"].values)
            exp2 += t[2] if len(t) > 2 else 0
            exp3 += t[3] if len(t) > 3 else 0
            exp4 += t[4] if len(t) > 4 else 0
        rows.append(dict(zip(KEY, combo),
                         genes=len(obs), p_rep=round(g["pass"].mean(), 4),
                         obs_ge2=int((obs >= 2).sum()), exp_ge2=round(exp2, 1),
                         obs_ge3=int((obs >= 3).sum()), exp_ge3=round(exp3, 2),
                         obs_4=int((obs >= 4).sum()), exp_4=round(exp4, 3)))
    r = pd.DataFrame(rows)
    tot = r[["genes", "obs_ge2", "exp_ge2", "obs_ge3", "exp_ge3", "obs_4", "exp_4"]].sum()
    print(f"\n==== {label} ====")
    print(r.to_string(index=False))
    print(f"  TOTAL  >=2 reps: observed {int(tot.obs_ge2)} vs expected by chance {tot.exp_ge2:.1f}"
          f" | >=3: {int(tot.obs_ge3)} vs {tot.exp_ge3:.1f} | all 4: {int(tot.obs_4)} vs {tot.exp_4:.2f}")
    r.to_csv(os.path.join(OUT, f"check_concordance_{label.split()[0]}.tsv"), sep="\t", index=False)

A = pd.read_csv(os.path.join(OUT, "mobile_transcripts_per_rep.tsv"), sep="\t")
run(A, "sum_log10BF", "snpsum (per-SNP sum)")
Bf = os.path.join(OUT, "mobile_gene_pooled_per_rep.tsv")
if os.path.exists(Bf):
    B = pd.read_csv(Bf, sep="\t")
    run(B, "log10BF", "pooled (gene-pooled)")
