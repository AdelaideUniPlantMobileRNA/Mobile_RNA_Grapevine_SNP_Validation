#!/usr/bin/env python3
# =====================================================================
# Check that every row's gene_id really contains the row's position (CHROM, POS)
# according to the GFF3. Used to detect the 09 row-order bug (gene ids assigned
# to the wrong rows).
# Usage: python check_gene_labels.py <gff3> <file.tsv> [more files...]
#   files need columns CHROM, POS, gene_id (header row); others are ignored.
# =====================================================================
import sys
import pandas as pd

gff = sys.argv[1]
genes = {}
with open(gff) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "gene":
            continue
        gid = [x[3:] for x in f[8].split(";") if x.startswith("ID=")]
        if gid:
            genes[gid[0]] = (f[0], int(f[3]), int(f[4]))
G = pd.DataFrame.from_dict(genes, orient="index", columns=["g_chrom", "g_start", "g_end"])
G.index.name = "gene_id"; G = G.reset_index()

for path in sys.argv[2:]:
    d = pd.read_csv(path, sep="\t", usecols=["CHROM", "POS", "gene_id"], dtype={"CHROM": str})
    d = d[d["gene_id"].notna()]
    if "source" in pd.read_csv(path, sep="\t", nrows=1).columns:   # per-sample 13 output
        d["source"] = pd.read_csv(path, sep="\t", usecols=["source"])["source"].loc[d.index]
    m = d.merge(G, on="gene_id", how="left")
    m["ok"] = (m["CHROM"] == m["g_chrom"]) & (m["POS"] >= m["g_start"]) & (m["POS"] <= m["g_end"])
    print(f"\n== {path}")
    print(f"   rows with gene_id: {len(m)}   label matches position: {int(m.ok.sum())} "
          f"({100 * m.ok.mean():.2f}%)")
    if "source" in m.columns:
        print(m.groupby("source")["ok"].agg(rows="size", correct="sum").to_string())
    bad = m[~m["ok"]].head(5)
    if len(bad):
        print("   examples of WRONG labels:")
        print(bad[["CHROM", "POS", "gene_id", "g_chrom", "g_start", "g_end"]].to_string(index=False))
