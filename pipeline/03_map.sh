#!/usr/bin/env bash
# Map each library in sample_sheet.csv to the four-variant reference, then
# reduce it to the compact per-sample files every later stage actually reads.
#
# --secondary=no is load-bearing: every downstream assay counts PRIMARY
# alignments only, and allowing secondary alignments would let one read be
# counted against two marker alleles at once.
#
# Reads are STREAMED from whatever the sheet's `fastq` column names -- a path,
# a glob, or zip:<archive>!<member>.  A 47 GB archive is never unpacked.
#
# KEEP_BAM=1 retains the whole-genome BAM (needed only by the genome-wide
# variant scan, stage 14).  Default is to keep the region BAM plus the depth
# tables, which is ~1% of the size and enough for stages 04-12.
set -euo pipefail
cd "$(dirname "$0")/.."
PIPE="pipeline"
REF="work/ref_probe.fa"
THREADS="${THREADS:-12}"
KEEP_BAM="${KEEP_BAM:-0}"
[ -s "$REF" ] || { echo "missing $REF -- run 02_build_ref.py first" >&2; exit 1; }
mkdir -p work/bam work/depth work/region

# Region set pulled once from config, so coordinates live in exactly one place.
# Written as a MERGED, SORTED bed and applied with -L: passing several regions
# to `samtools view` instead would emit reads in region order (not coordinate
# order, so the result cannot be indexed) and would emit a read TWICE when two
# requested regions overlap, silently inflating every count downstream.
BED="work/regions.bed"
python - "$BED" <<'EOF'
import sys
sys.path.insert(0, 'pipeline'); import config as C
PAD = 5000
iv = []
for lo, hi in (C.CASSETTE_WINDOW, C.PROPHAGE, C.TFA,
               C.GENES['lamB']['span'], C.GENES['malT']['span']):
    iv.append((max(0, lo - PAD), hi + PAD))
iv.sort()
merged = [list(iv[0])]
for lo, hi in iv[1:]:
    if lo <= merged[-1][1]:
        merged[-1][1] = max(merged[-1][1], hi)
    else:
        merged.append([lo, hi])
lines = ['%s\t%d\t%d' % (C.CONTIG, lo, hi) for lo, hi in merged]
# whole decoy contigs: short, and 04 needs every read that landed on them
import json, os
meta = json.load(open(os.path.join(C.WORK, 'ref_probe.json')))
for a in C.ALLELES:
    lines.append('%s\t0\t%d' % (a, meta[a]['length']))
open(sys.argv[1], 'w').write('\n'.join(lines) + '\n')
print('regions:', len(lines), 'bp on %s:' % C.CONTIG,
      sum(h - l for l, h in merged))
EOF

for S in $(python -c "import sys; sys.path.insert(0,'pipeline'); import config as C; print(' '.join(C.SAMPLES))"); do
  if [ -s "work/depth/${S}.chrom1kb.txt" ] && [ -s "work/region/${S}.bam.bai" ]; then
    echo "skip $S (done)"; continue
  fi
  echo "=== $S ==="
  FQCMD=$(python -c "import sys; sys.path.insert(0,'pipeline'); import config as C; print(C.fastq_cmd('$S'))")
  eval "$FQCMD" \
    | minimap2 -ax map-ont -t "$THREADS" --secondary=no -K 200M "$REF" - 2>"work/depth/${S}.mm2.log" \
    | samtools sort -@ 3 -m 1G -o "work/bam/${S}.bam" -
  samtools index -@ 4 "work/bam/${S}.bam"

  # chromosome depth in 1 kb bins -- copy number (06) and the deletion scan
  samtools depth -a -Q 0 -r "$(python -c "import sys; sys.path.insert(0,'pipeline'); import config as C; print(C.CONTIG)")" "work/bam/${S}.bam" \
    | awk 'BEGIN{OFS="\t"} {b=int(($2-1)/1000); s[b]+=$3; n[b]++}
           END{for(i in s) print i*1000+1, s[i]/n[i]}' \
    | sort -k1,1n > "work/depth/${S}.chrom1kb.txt"
  samtools flagstat "work/bam/${S}.bam" > "work/depth/${S}.flagstat.txt"

  # compact region BAM: everything stages 04-12 look at
  samtools view -b -L "$BED" -o "work/region/${S}.bam" "work/bam/${S}.bam"
  samtools index "work/region/${S}.bam"

  if [ "$KEEP_BAM" = "0" ]; then rm -f "work/bam/${S}.bam" "work/bam/${S}.bam.bai"; fi
  echo "done $S"
done
echo "mapping complete"
