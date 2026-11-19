#!/usr/bin/env python3
# =====================================================================
# Real-data spot check of the Bayesian step (13 v3), independent of our code.
#
# For randomly chosen SNP markers of one heterograft sample:
#  1) re-count reads DIRECTLY from the BAM files with samtools mpileup
#     (a different program from the bcftools mpileup used in the pipeline):
#       heterograft N, n ; local homograft Nh1, nh1 ; distal homograft Nh2, nh2
#  2) re-compute the Bayes factor with baymobil's OFFICIAL single-value function
#     run_bayes(Nh1, nh1, Nh2, nh2, N, n, nmax) from the published package
#  3) compare both with what the pipeline wrote in per_sample/<sample>.snp.tsv
#  4) check the donor allele is really donor-private in the marker databases
#  5) check the gene sum in mobile_transcripts_per_rep.tsv = sum of the SNP values
#
# Run from the project directory:
#   python scripts/19_spotcheck.py [sample=SRL1] [n_markers=20] [seed=1]
# =====================================================================
import sys, os, re, subprocess, random
import numpy as np
import pandas as pd

SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "SRL1"
NMARK  = int(sys.argv[2]) if len(sys.argv) > 2 else 20
SEED   = int(sys.argv[3]) if len(sys.argv) > 3 else 1
REF  = "/scratchdata1/users/a1709557/Own_RNA/Reference_genomes/Grape_reference_assemblies/T2T_ref.fasta"
OUT  = "results/mobile_v3"
CAP, NMAX = 2000, 10
BASES = "ACGT"

from baymobil.baymobil import run_bayes          # official published entry point

het = pd.read_csv("metadata/hetero_samples.tsv", sep="\t")
hom = pd.read_csv("metadata/step2_samples.tsv", sep="\t")
s = het[het["sample"] == SAMPLE].iloc[0]
local_bams  = [f"results/star/{x}.bam" for x in
               hom[(hom.group == s.local_geno) & (hom.tissue == s.tissue)]["sample"]]
distal_bams = [f"results/star/{x}.bam" for x in hom[hom.group == s.mobile_geno]["sample"]]
het_bam = f"results/star_hetero/{SAMPLE}.bam"

# ---- pick markers: half with log10BF >= 1, half random (fixed seed = reproducible) ----
d = pd.read_csv(f"{OUT}/per_sample/{SAMPLE}.snp.tsv", sep="\t", dtype={"CHROM": str})
pos_ = d[d.log10BF >= 1]; rest = d[d.log10BF < 1]
k = min(len(pos_), NMARK // 2)
pick = pd.concat([pos_.sample(k, random_state=SEED),
                  rest.sample(NMARK - k, random_state=SEED)]).reset_index(drop=True)

def parse_bases(ref, bases):
    """count A,C,G,T in one samtools mpileup base column."""
    c = dict.fromkeys(BASES, 0)
    i = 0
    while i < len(bases):
        ch = bases[i]
        if ch == "^":                      # read start + mapping-quality char
            i += 2; continue
        if ch == "$":
            i += 1; continue
        if ch in "+-":                     # indel: skip its sequence
            m = re.match(r"[+-](\d+)", bases[i:])
            i += len(m.group(0)) + int(m.group(1)); continue
        if ch in ".,":
            c[ref.upper()] += 1
        elif ch.upper() in BASES:
            c[ch.upper()] += 1
        i += 1                             # '*', '#', 'N' etc. are not counted
    return c

def samtools_counts(chrom, pos, bams):
    cmd = ["samtools", "mpileup", "-f", REF, "-r", f"{chrom}:{pos}-{pos}",
           "-q", "20", "-Q", "20", "-B", "-d", "10000"] + bams
    line = subprocess.run(cmd, capture_output=True, text=True).stdout.strip()
    tot = dict.fromkeys(BASES, 0)
    if not line:
        return tot
    f = line.split("\t"); ref = f[2]
    for j in range(len(bams)):
        b = f[4 + 3 * j] if 4 + 3 * j < len(f) else ""
        for base, v in parse_bases(ref, b).items():
            tot[base] += v
    return tot

def cap(Nv, nv):
    if Nv > CAP:
        return CAP, int(np.round(nv * CAP / Nv))
    return Nv, nv

diag = pd.read_csv("results/dipcall/diagnostic_sites.tsv", sep="\t", dtype={"CHROM": str},
                   usecols=["CHROM", "POS", "ruggeri_private", "shiraz_private"])
pcol = "ruggeri_private" if s.direction == "ruggeri" else "shiraz_private"
mdb = pd.read_csv("results/dipcall/master_diagnostic_db_v3.tsv", sep="\t", dtype={"CHROM": str},
                  usecols=["CHROM", "POS", "tissue", "direction", "allele"])

rows = []
for _, r in pick.iterrows():
    a = r.allele
    h = samtools_counts(r.CHROM, r.POS, [het_bam])
    l = samtools_counts(r.CHROM, r.POS, local_bams)
    t = samtools_counts(r.CHROM, r.POS, distal_bams)
    N, n = sum(h.values()), h[a]
    Nh1, nh1 = cap(sum(l.values()), l[a])
    Nh2, nh2 = cap(sum(t.values()), sum(t.values()) - t[a])
    bf_pipeline_inputs = run_bayes(int(r.Nh1), int(r.nh1), int(r.Nh2), int(r.nh2),
                                   int(r.N), int(r.n), NMAX)[2]
    bf_independent = run_bayes(int(Nh1), int(nh1), int(Nh2), int(nh2), int(N), int(n), NMAX)[2] \
        if N > 0 else np.nan
    if r.source == "genome":
        g = diag[(diag.CHROM == r.CHROM) & (diag.POS == r.POS)][pcol].astype(str)
        ok_allele = any(a in x.split(",") for x in g)
    else:
        ok_allele = len(mdb[(mdb.CHROM == r.CHROM) & (mdb.POS == r.POS) & (mdb.tissue == s.tissue) &
                            (mdb.direction == s.direction) & (mdb.allele == a)]) > 0
    rows.append(dict(CHROM=r.CHROM, POS=r.POS, allele=a, source=r.source,
                     N_pipe=r.N, N_sam=N, n_pipe=r.n, n_sam=n,
                     Nh1_pipe=r.Nh1, Nh1_sam=Nh1, nh1_pipe=r.nh1, nh1_sam=nh1,
                     Nh2_pipe=r.Nh2, Nh2_sam=Nh2, nh2_pipe=r.nh2, nh2_sam=nh2,
                     BF_pipe=round(r.log10BF, 3), BF_official=round(bf_pipeline_inputs, 3),
                     BF_from_samtools=round(bf_independent, 3), allele_in_db=ok_allele))
    print(f"   checked {r.CHROM}:{r.POS}", flush=True)

T = pd.DataFrame(rows)
T.to_csv(f"{OUT}/spotcheck_{SAMPLE}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
print(f"\n==== spot check: {SAMPLE} ({s.graft} {s.tissue}, {s.direction}), {len(T)} markers ====")
print(T.to_string(index=False))

def close(x, y, tol):
    return (abs(x - y) <= tol * np.maximum(1, np.abs(y)))
print("\n==== summary (pipeline vs independent) ====")
print(f"  heterograft total reads N  equal within 2%: {close(T.N_pipe, T.N_sam, 0.02).sum()}/{len(T)}")
print(f"  donor-allele reads n       equal within 2%: {close(T.n_pipe, T.n_sam, 0.02).sum()}/{len(T)}")
print(f"  local homograft Nh1        equal within 2%: {close(T.Nh1_pipe, T.Nh1_sam, 0.02).sum()}/{len(T)}")
print(f"  local homograft nh1        equal within 2%: {close(T.nh1_pipe, T.nh1_sam, 0.02).sum()}/{len(T)}")
print(f"  distal homograft Nh2       equal within 2%: {close(T.Nh2_pipe, T.Nh2_sam, 0.02).sum()}/{len(T)}")
print(f"  pipeline BF == official run_bayes on same inputs (|diff|<0.001): "
      f"{(abs(T.BF_pipe - T.BF_official) < 1e-3).sum()}/{len(T)}")
print(f"  same call (log10BF>=1) from samtools counts: "
      f"{((T.BF_pipe >= 1) == (T.BF_from_samtools >= 1)).sum()}/{len(T)}")
print(f"  donor allele found in marker database: {T.allele_in_db.sum()}/{len(T)}")

# ---- gene sums ----
pr = pd.read_csv(f"{OUT}/mobile_transcripts_per_rep.tsv", sep="\t")
pr = pr[pr["sample"] == SAMPLE].set_index("gene_id")["sum_log10BF"]
genes = pick.gene_id.unique()
own = d[d.gene_id.isin(genes)].groupby("gene_id")["log10BF"].sum()
ok = sum(abs(own[g] - pr.get(g, np.nan)) < 1e-6 for g in genes)
print(f"  gene sum in mobile_transcripts_per_rep == sum of its SNP values: {ok}/{len(genes)}")
print(f"\n  -> {OUT}/spotcheck_{SAMPLE}.tsv")
