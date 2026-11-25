#!/usr/bin/env bash
# Full pipeline.  Edit pipeline/config_inputs.py (where the data is) and
# pipeline/config.py (the biological system); nothing else needs changing.
#
#   ./pipeline/run_all.sh            # everything
#   ./pipeline/run_all.sh 07         # from stage 07 onward
#   KEEP_BAM=1 ./pipeline/run_all.sh # keep whole-genome BAMs (needed by 14)
set -euo pipefail
cd "$(dirname "$0")/.."
FROM="${1:-00}"
run() {
  local n="${1%%_*}"; [ "$n" \< "$FROM" ] && { echo "skip $1"; return; }
  echo "=== $1 ==="
  case "$1" in *.sh) bash "pipeline/$1";; *.py) python "pipeline/$1";; esac
}
run 00_setup_env.sh          # once per machine; comment out afterwards
run 00_sample_sheet.py       # plate record + archives -> sample_sheet.csv
run 02_build_ref.py          # masked backbone + allele decoys (mapping reference)
run 02b_build_locus_ref.py   # four-variant locus reference (identical anchors)
run 03_map.sh                # map, then reduce to region BAM + 1 kb depth
run 01_qc.py                 # read_stats.csv
run 04_locus_assign.py       # marker allele fractions, direct read assignment
run 05_lysis_junction.py     # cassette-independent lysis status + unmixing
run 06_copy_number.py        # prophage copies per chromosome
run 07_sv_scan.py            # genome-wide structural events
run 08_extract_inserts.py
run 09_receptor_alleles.py   # lamB lesion classes
run 10_is_junctions.py       # IS family per junction
run 11_activator.py          # malT
run 12_tfa_recomb.py         # prophage tail-fibre recombinants
run 13_parental_catalogue.sh # sites that differ between the ancestors
run 14_variant_scan.sh       # genome-wide joint call (needs KEEP_BAM=1)
run 15_variant_filter.py     # parental / differential / repeat filters
run 30_tables.py             # result CSVs
echo "pipeline complete"
