#!/usr/bin/env bash
# Catalogue of sites where the OTHER parental assemblies differ from the
# mapping reference. These are not evolved mutations -- they are the starting
# differences between the competing strains, and they must be subtracted from
# the genome-wide scan or every mixed flask looks hypermutated.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p work/par
eval "$(python -c "
import sys; sys.path.insert(0,'pipeline'); import config as C
print('REFFA=%r; CONTIG=%r; ANC=%r' % (C.REF_FASTA, C.CONTIG, C.ANCESTOR_DIR))
print('PARENTS=(%s)' % ' '.join(repr(p) for p in C.PARENTS[1:]))")"
# seqkit, not `samtools faidx`: the reference often lives on a read-only
# mount and faidx insists on writing a .fai next to it.
seqkit grep -n -p "$CONTIG" "$REFFA" > work/par/ref_contig.fa
for P in "${PARENTS[@]}"; do
  Q="$ANC/$P/${P}_1_reference.fasta"
  [ -s "$Q" ] || { echo "missing $Q" >&2; continue; }
  nucmer --maxmatch -p "work/par/${P}_" work/par/ref_contig.fa "$Q" 2>/dev/null
  delta-filter -1 "work/par/${P}_.delta" > "work/par/${P}_.1delta"
  show-snps -Clr "work/par/${P}_.1delta" > "work/par/${P}_.snps"
done
awk 'NR>4 && $2!="." && $3!="." {print $1}' work/par/*.snps | sort -un > work/par/parental_sites.txt
echo "parental SNV sites: $(wc -l < work/par/parental_sites.txt)"
