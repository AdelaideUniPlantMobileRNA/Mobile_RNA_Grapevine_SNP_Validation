#!/usr/bin/env python3
# =====================================================================
# Step 4b (v3): combine per-sample results and report TWO summaries side by side.
#   A  per-SNP sum (baymobil original, Tomkins et al. 2022):
#        gene log10BF in a replicate = SUM of its SNP log10BF   (<sample>.snp.tsv)
#   B  gene-pooled:
#        gene log10BF in a replicate = ONE Bayes factor on counts pooled over SNPs
#        >= 300 bp apart                                       (<sample>.gene.tsv)
#   For both: replicate passes if gene log10BF >= 1; gene is MOBILE if it passes
#   in >= MIN_REPS replicates of (graft, tissue, direction).
# Outputs in <out_dir>:
#   mobile_transcripts.tsv         A (same format as before; used by 16_inspect_mobile.py)
#   mobile_transcripts_pooled.tsv  B
#   mobile_transcripts_both.tsv    A and B side by side, one row per gene x combination
# Usage: python 13b_aggregate_mobile_v3.py <out_dir> [min_reps=2] [exclude=S1,S2]
# =====================================================================
import sys, os, glob
import pandas as pd

OUT = sys.argv[1]
MIN_REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 2
EXCLUDE = set(x for x in (sys.argv[3].split(",") if len(sys.argv) > 3 else []) if x)
THR = 1.0
KEY = ["graft", "tissue", "direction", "gene_id"]

def load(suffix):
    fs = sorted(glob.glob(os.path.join(OUT, "per_sample", f"*.{suffix}.tsv")))
    fs = [f for f in fs if os.path.basename(f).replace(f".{suffix}.tsv", "") not in EXCLUDE]
    print(f">> {suffix}: using {len(fs)} sample files (excluded: {sorted(EXCLUDE) or 'none'})")
    return pd.concat([pd.read_csv(f, sep="\t", dtype={"CHROM": str}) for f in fs], ignore_index=True)

# ---- A: per-SNP sum ----
snp = load("snp")
snp.to_csv(os.path.join(OUT, "mobile_snp_log10BF.tsv"), sep="\t", index=False)
txA = (snp.groupby(KEY[:3] + ["rep", "sample", "gene_id"], as_index=False)
          .agg(n_snp=("log10BF", "size"), sum_log10BF=("log10BF", "sum"),
               sum_N2max=("N2max", "sum")))
txA["rep_mobile"] = txA["sum_log10BF"] >= THR
txA.to_csv(os.path.join(OUT, "mobile_transcripts_per_rep.tsv"), sep="\t", index=False)
A = (txA.groupby(KEY, as_index=False)
        .agg(n_reps=("rep", "nunique"), n_reps_mobile=("rep_mobile", "sum"),
             median_n_snp=("n_snp", "median"), max_sum_log10BF=("sum_log10BF", "max")))
A["mobile"] = A["n_reps_mobile"] >= MIN_REPS
A.sort_values(["mobile", "n_reps_mobile", "max_sum_log10BF"], ascending=False) \
 .to_csv(os.path.join(OUT, "mobile_transcripts.tsv"), sep="\t", index=False)

# ---- B: gene-pooled ----
gen = load("gene")
gen["rep_mobile"] = gen["log10BF"] >= THR
gen.to_csv(os.path.join(OUT, "mobile_gene_pooled_per_rep.tsv"), sep="\t", index=False)
B = (gen.groupby(KEY, as_index=False)
        .agg(n_reps=("rep", "nunique"), n_reps_mobile=("rep_mobile", "sum"),
             median_snp_used=("n_snp_used", "median"), max_log10BF=("log10BF", "max")))
B["mobile"] = B["n_reps_mobile"] >= MIN_REPS
B.sort_values(["mobile", "n_reps_mobile", "max_log10BF"], ascending=False) \
 .to_csv(os.path.join(OUT, "mobile_transcripts_pooled.tsv"), sep="\t", index=False)

# ---- side by side ----
both = A.rename(columns={"n_reps_mobile": "reps_pass_snpsum", "max_sum_log10BF": "max_log10BF_snpsum",
                         "mobile": "mobile_snpsum"}) \
        .merge(B[KEY + ["n_reps_mobile", "median_snp_used", "max_log10BF", "mobile"]]
               .rename(columns={"n_reps_mobile": "reps_pass_pooled", "max_log10BF": "max_log10BF_pooled",
                                "mobile": "mobile_pooled"}), on=KEY, how="outer")
for c in ["mobile_snpsum", "mobile_pooled"]:
    both[c] = both[c].fillna(False).astype(bool)
both["called_by"] = both.apply(lambda r: "both" if r.mobile_snpsum and r.mobile_pooled else
                               ("snpsum_only" if r.mobile_snpsum else
                                ("pooled_only" if r.mobile_pooled else "")), axis=1)
both = both.sort_values(["called_by", "graft", "tissue", "direction", "gene_id"])
both.to_csv(os.path.join(OUT, "mobile_transcripts_both.tsv"), sep="\t", index=False)

pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
print(f"\n==== MOBILE (>= {MIN_REPS} reps with gene log10BF >= {THR:g}) : two summaries ====")
print("  gene x combination tested:", len(both))
tab = both.groupby(["graft", "tissue", "direction"]).agg(
    snpsum=("mobile_snpsum", "sum"), pooled=("mobile_pooled", "sum"),
    both=("called_by", lambda s: (s == "both").sum())).reset_index()
print(tab.to_string(index=False))
print("  TOTAL  per-SNP sum:", int(both.mobile_snpsum.sum()),
      "  gene-pooled:", int(both.mobile_pooled.sum()),
      "  called by both:", int((both.called_by == "both").sum()))
print("\n  called genes:")
print(both[both.called_by != ""][KEY + ["called_by", "reps_pass_snpsum", "reps_pass_pooled",
      "max_log10BF_snpsum", "max_log10BF_pooled"]].to_string(index=False))
print("  ->", os.path.join(OUT, "mobile_transcripts_both.tsv"))
