#!/usr/bin/env bash
# Joint genome-wide variant scan across all samples, then subtraction of the
# parental catalogue from 13.
#
# The per-file depth cap matters enormously: at -d 2000 this does not finish on
# 11 ONT BAMs (killed after 2 h with no output). MPILEUP_D in config.py is the
# knob; 150 is ample for allele-frequency estimation at these depths.
set -euo pipefail
cd "$(dirname "$0")/.."
eval "$(python -c "
import sys; sys.path.insert(0,'pipeline'); import config as C
print('REFFA=%r; D=%d; Q=%d; DP=%d' % (C.MAP_FASTA, C.MPILEUP_D, C.MIN_QUAL, C.MIN_DP))
print('BAMS=(%s)' % ' '.join(repr(C.bam(s)) for s in C.SAMPLES))")"
bcftools mpileup -f "$REFFA" -d "$D" -q 20 -Q 10 -L "$D" -B \
    -a FORMAT/AD,FORMAT/DP -Ou "${BAMS[@]}" 2> work/par/mpileup.log \
  | bcftools call -m -v -Ou 2>> work/par/mpileup.log \
  | bcftools filter -i "QUAL>=$Q && INFO/DP>=$DP" -Oz -o work/par/all.vcf.gz
bcftools index -f work/par/all.vcf.gz
# drop sites already different between the parents
awk '{print "'"$(python -c "
import sys; sys.path.insert(0,'pipeline'); import config as C; print(C.CONTIG)")"'\t"$1}' \
    work/par/parental_sites.txt > work/par/parental.tab
bcftools view -T ^work/par/parental.tab -Oz -o work/par/evolved.vcf.gz work/par/all.vcf.gz
bcftools index -f work/par/evolved.vcf.gz
echo "all sites:      $(bcftools view -H work/par/all.vcf.gz | wc -l)"
echo "after parental: $(bcftools view -H work/par/evolved.vcf.gz | wc -l)"
