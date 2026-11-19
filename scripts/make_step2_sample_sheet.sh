#!/bin/bash
# Build metadata/step2_samples.tsv from the trimmed SS and RR FASTQ folders.
# Sample name convention: <group><tissue><rep>, e.g. SSL2 = SS, tissue L, rep 2.
set -euo pipefail

SS_DIR=/scratchdata1/users/a1709557/Own_RNA/02.TrimmedData/T1T2/SS_Samples
RR_DIR=/scratchdata1/users/a1709557/Own_RNA/02.TrimmedData/T1T2/RR_Samples
OUT=metadata/step2_samples.tsv

mkdir -p metadata
printf 'sample\tgroup\ttissue\trep\tR1\tR2\n' > "$OUT"
for d in "$SS_DIR" "$RR_DIR"; do
  for r1 in "$d"/*.R1.trimmed.fq.gz; do
    [[ -e "$r1" ]] || continue
    s=$(basename "$r1" .R1.trimmed.fq.gz)
    r2=${r1/.R1.trimmed.fq.gz/.R2.trimmed.fq.gz}
    grp=${s:0:2}                      # SS or RR
    rest=${s:2}                       # e.g. L2, PY4, B1, R3
    tissue=$(sed 's/[0-9]*$//' <<<"$rest")
    rep=$(grep -o '[0-9]*$' <<<"$rest")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$grp" "$tissue" "$rep" "$r1" "$r2" >> "$OUT"
  done
done
echo ">> wrote $OUT  ($(($(wc -l < "$OUT")-1)) samples)"
column -t "$OUT"
