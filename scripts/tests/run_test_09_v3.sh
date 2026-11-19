#!/bin/bash
# Self-contained verification of 09_build_master_db_v3.R (L1-L4 + strand-aware L3).
# v3 extra case: pos 700, Ruggeri allele expressed ONLY in Ruggeri ROOT (RRR),
# (A->C, not an editing change); SSL covered without it -> must become a ruggeri
# marker for tissue L with donor_tissues = R.
# Toy count tables exercise: a clean genome+RNA marker, an A->G RNA-editing site
# (removed), a low-reproducibility call (removed), a strand-biased call, a
# recipient-also-expressed allele, and a REAL T->C SNP on a (+)-strand gene
# (must be KEPT, because editing on a + gene is A->G, not T->C).
#   bash run_test_09.sh        (needs Rscript with data.table)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
R9="$HERE/../09_build_master_db_v3.R"
[[ -f "$R9" ]] || { echo "MISSING: $R9"; exit 1; }
WD="$HERE/_t09v3"; rm -rf "$WD"; mkdir -p "$WD/cnt"

cat > "$WD/step2_samples.tsv" <<'EOF'
sample	group	tissue	rep
SSL1	SS	L	1
SSL2	SS	L	2
RRL1	RR	L	1
RRL2	RR	L	2
RRR1	RR	R	1
RRR2	RR	R	2
EOF

cat > "$WD/diagnostic_sites.tsv" <<'EOF'
CHROM	POS	REF	ALT	shiraz_alleles	ruggeri_alleles	shiraz_private	ruggeri_private	detects_shiraz_in_rootstock	detects_ruggeri_in_scion	is_indel	gene_id	gene_name	is_exonic
c1	100	A	C	C	A	C	A	TRUE	TRUE	FALSE	g1	GENE1	TRUE
EOF

# one (+)-strand gene spanning all toy positions
printf 'c1\tx\tgene\t1\t650\t.\t+\t.\tID=g1;Name=GENE1\nc1\tx\tgene\t651\t10000\t.\t+\t.\tID=g2;Name=GENE2\n' > "$WD/toy.gff3"

cat > "$WD/cnt/SSL1.cnt.tsv" <<'EOF'
c1	100	A	C,<*>	0,20,0	0,10,0	0,10,0
c1	200	A	G,<*>	15,5,0	8,3,0	7,2,0
c1	300	C	T,<*>	10,10,0	5,5,0	5,5,0
c1	400	G	A,<*>	0,12,0	0,12,0	0,0,0
c1	500	A	C,<*>	0,20,0	0,10,0	0,10,0
c1	600	T	C,<*>	0,20,0	0,10,0	0,10,0
c1	700	A	<*>	40,0	20,0	20,0
EOF
cat > "$WD/cnt/SSL2.cnt.tsv" <<'EOF'
c1	100	A	C,<*>	0,18,0	0,9,0	0,9,0
c1	200	A	G,<*>	14,4,0	7,2,0	7,2,0
c1	300	C	<*>	25,0	13,0	12,0
c1	400	G	A,<*>	0,11,0	0,11,0	0,0,0
c1	500	A	C,<*>	0,19,0	0,9,0	0,10,0
c1	600	T	C,<*>	0,18,0	0,9,0	0,9,0
c1	700	A	<*>	35,0	18,0	17,0
EOF
cat > "$WD/cnt/RRL1.cnt.tsv" <<'EOF'
c1	100	A	<*>	30,0	15,0	15,0
c1	200	A	<*>	26,0	13,0	13,0
c1	300	C	<*>	22,0	11,0	11,0
c1	400	G	<*>	20,0	10,0	10,0
c1	500	A	C,<*>	15,15,0	8,8,0	7,7,0
c1	600	T	<*>	30,0	15,0	15,0
EOF
cat > "$WD/cnt/RRL2.cnt.tsv" <<'EOF'
c1	100	A	<*>	28,0	14,0	14,0
c1	200	A	<*>	24,0	12,0	12,0
c1	300	C	<*>	20,0	10,0	10,0
c1	400	G	<*>	18,0	9,0	9,0
c1	500	A	C,<*>	14,16,0	7,8,0	7,8,0
c1	600	T	<*>	28,0	14,0	14,0
EOF

cat > "$WD/cnt/RRR1.cnt.tsv" <<'EOF'
c1	700	A	C,<*>	0,30,0	0,15,0	0,15,0
EOF
cat > "$WD/cnt/RRR2.cnt.tsv" <<'EOF'
c1	700	A	C,<*>	0,25,0	0,12,0	0,13,0
EOF

Rscript "$R9" "$WD" "$WD/cnt" "$WD/step2_samples.tsv" "$WD/toy.gff3" 10 3 0.10 2 10 0.02
OUT="$WD/master_diagnostic_db_v3.tsv"
echo "----- produced markers -----"; cut -f1,2,4,5,6,10,12 "$OUT"

has () { awk -F'\t' -v p="$1" -v d="$2" -v a="$3" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
    $h["POS"]==p&&$h["direction"]==d&&$h["allele"]==a{print $h["tier"]}' "$OUT"; }
rows_at () { awk -F'\t' -v p="$1" 'NR>1 && $2==p{c++} END{print c+0}' "$OUT"; }

F=0
[[ "$(has 100 shiraz  C)" == "genome+RNA" ]] || { echo "FAIL 100 shiraz C";  F=1; }
[[ "$(has 100 ruggeri A)" == "genome+RNA" ]] || { echo "FAIL 100 ruggeri A"; F=1; }
[[ "$(has 400 ruggeri G)" == "RNA-only"   ]] || { echo "FAIL 400 ruggeri G"; F=1; }
[[ "$(has 500 ruggeri A)" == "RNA-only"   ]] || { echo "FAIL 500 ruggeri A"; F=1; }
[[ "$(has 600 shiraz  C)" == "RNA-only"   ]] || { echo "FAIL 600 shiraz C (real T->C on + gene wrongly removed)"; F=1; }
[[ "$(rows_at 200)" == 0 ]] || { echo "FAIL: A->G editing site 200 leaked"; F=1; }
[[ "$(rows_at 300)" == 0 ]] || { echo "FAIL: low-repro 300 leaked"; F=1; }
[[ -z "$(has 400 shiraz A)" ]] || { echo "FAIL: strand-biased A called shiraz-private"; F=1; }
[[ -z "$(has 500 shiraz C)" ]] || { echo "FAIL: recipient-expressed C called shiraz-private"; F=1; }
T7=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["POS"]==700&&$h["tissue"]=="L"&&$h["direction"]=="ruggeri"&&$h["allele"]=="C"{print $h["donor_tissues"]}' "$OUT")
[[ "$T7" == "R" ]] || { echo "FAIL 700: root-only Ruggeri allele not kept for leaf (got '$T7')"; F=1; }
# gene labels must match position (g1 = 1-650, g2 = 651-10000), for every row
BAD=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
   { want=($h["POS"]<=650)?"g1":"g2"; if($h["gene_id"]!=want) n++ } END{print n+0}' "$OUT")
[[ "$BAD" == 0 ]] || { echo "FAIL: $BAD rows carry the wrong gene_id (row-order bug)"; F=1; }
if [[ "$F" == 0 ]]; then echo "========= TEST PASSED ========="; else echo "========= TEST FAILED ========="; exit 1; fi
