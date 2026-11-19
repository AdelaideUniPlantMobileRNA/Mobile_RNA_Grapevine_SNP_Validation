#!/bin/bash
# Build metadata/hetero_samples.tsv for the heterograft samples.
# Each folder is one graft x sampled-tissue-genotype; the local genotype, the
# mobile (distal) genotype we look for, and the direction are fixed per folder.
# Reads are aligned to T2T (common backbone), same as the homografts.
#   local_geno = the genotype the tissue belongs to (its tissue-matched homograft
#                is the background q reference)
#   mobile_geno / direction = the OTHER genotype's private alleles to detect
set -euo pipefail

BASE=/scratchdata1/users/a1709557/Own_RNA/02.TrimmedData/T1T2
OUT=metadata/hetero_samples.tsv
mkdir -p metadata
printf 'sample\tgraft\tlocal_geno\tmobile_geno\tdirection\ttissue\trep\tR1\tR2\n' > "$OUT"

add () {  # folder  graft(SR/RS)  local_geno(SS/RR)  mobile_geno(SS/RR)  direction(shiraz/ruggeri)
  local folder=$1 graft=$2 lg=$3 mg=$4 dir=$5
  for r1 in "$BASE/$folder"/*.R1.trimmed.fq.gz; do
    [[ -e "$r1" ]] || continue
    s=$(basename "$r1" .R1.trimmed.fq.gz)
    r2=${r1/.R1.trimmed.fq.gz/.R2.trimmed.fq.gz}
    rest=${s#$graft}                                   # strip SR/RS -> tissue+rep (e.g. B1, PY3, R1, R12)
    tissue=$(sed 's/[0-9]*$//' <<<"$rest")
    rep=$(grep -o '[0-9]*$' <<<"$rest")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$graft" "$lg" "$mg" "$dir" "$tissue" "$rep" "$r1" "$r2" >> "$OUT"
  done
}

#    folder         graft local mobile direction
add  SRscion         SR    SS    RR     ruggeri    # Shiraz scion tissue  -> detect Ruggeri (up)
add  SRR_Samples     SR    RR    SS     shiraz     # Ruggeri root  tissue -> detect Shiraz  (down)
add  RSscion         RS    RR    SS     shiraz     # Ruggeri scion tissue -> detect Shiraz  (up)
add  RSR_Samples     RS    SS    RR     ruggeri    # Shiraz root   tissue -> detect Ruggeri (down)

echo ">> wrote $OUT  ($(($(wc -l < "$OUT")-1)) heterograft samples)"
column -t "$OUT" | cut -c1-120
