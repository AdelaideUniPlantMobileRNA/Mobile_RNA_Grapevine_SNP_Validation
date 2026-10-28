#!/bin/bash
# Self-contained verification of 02_build_diagnostic_sites.R.
# Writes its own toy data + truth table, runs the script, diffs the result.
# Usage:  bash run_test.sh      (needs Rscript with data.table on PATH)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
R2="$HERE/../02_build_diagnostic_sites.R"
[[ -f "$R2" ]] || { echo "MISSING: $R2  (put 02_build_diagnostic_sites.R in scripts/)"; exit 1; }

# ---- toy merged_genotypes.tsv : CHROM POS REF ALT shiraz_GT ruggeri_GT ----
cat > "$HERE/merged_genotypes.tsv" <<'EOF'
c1	150	A	C	1|1	0|0
c1	160	A	T	1|.	0|0
c1	420	G	T	0|1	1|1
c1	450	G	A,T	1|1	2|2
c1	600	C	G	./.	1|1
c1	720	A	G	1|1	1|1
c1	730	A	G	0|1	0|1
c1	800	A	T	1|1	./.
c1	1500	A	T	1|1	./.
EOF

# ---- primary (hap1) coverage: both genotypes cover c1:1-1000 (0-based half-open) ----
printf 'c1\t0\t1000\n' > "$HERE/shiraz.hap1.bed"
printf 'c1\t0\t1000\n' > "$HERE/ruggeri.hap1.bed"

# ---- toy gff3 : g1 100-500 (exons 100-200,400-500) ; g2 700-900 (exon 700-750) ----
cat > "$HERE/toy.gff3" <<'EOF'
c1	x	gene	100	500	.	+	.	ID=g1;Name=GENE1
c1	x	exon	100	200	.	+	.	ID=g1.e1;Parent=g1
c1	x	exon	400	500	.	+	.	ID=g1.e2;Parent=g1
c1	x	gene	700	900	.	+	.	ID=g2;Name=GENE2
c1	x	exon	700	750	.	+	.	ID=g2.e1;Parent=g2
EOF

# ---- hand-checked truth table (shipping output schema) ----
cat > "$HERE/expected_diagnostic_sites.tsv" <<'EOF'
CHROM	POS	REF	ALT	shiraz_alleles	ruggeri_alleles	shiraz_private	ruggeri_private	detects_shiraz_in_rootstock	detects_ruggeri_in_scion	is_indel	gene_id	gene_name	is_exonic
c1	150	A	C	C	A	C	A	TRUE	TRUE	FALSE	g1	GENE1	TRUE
c1	160	A	T	T	A	T	A	TRUE	TRUE	FALSE	g1	GENE1	TRUE
c1	420	G	T	G,T	T	G	NA	TRUE	FALSE	FALSE	g1	GENE1	TRUE
c1	450	G	A,T	A	T	A	T	TRUE	TRUE	FALSE	g1	GENE1	TRUE
c1	600	C	G	C	G	C	G	TRUE	TRUE	FALSE	NA	NA	FALSE
c1	800	A	T	T	A	T	A	TRUE	TRUE	FALSE	g2	GENE2	FALSE
EOF

echo ">> running 02 on toy data ..."
Rscript "$R2" "$HERE" "$HERE/toy.gff3"

echo "----- diff (expected vs produced); no lines between the markers = PASS -----"
if diff <(sort "$HERE/expected_diagnostic_sites.tsv") <(sort "$HERE/diagnostic_sites.tsv"); then
  echo "========================  TEST PASSED  ========================"
else
  echo "========================  TEST FAILED  ========================"; exit 1
fi
