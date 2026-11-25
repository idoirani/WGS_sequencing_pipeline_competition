#!/usr/bin/env python
"""Per-sample read statistics -> read_stats.csv.

Run AFTER 03_map.sh: pct_mapped and depth come from the mapping outputs.
Length/quality statistics are computed by streaming the reads through seqkit,
so a library inside a zip is never unpacked.
"""
import os, sys, csv, subprocess, statistics
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

rows = []
for s in C.SAMPLES:
    st_out = subprocess.run('%s | seqkit stats -a -T -' % C.fastq_cmd(s),
                            shell=True, capture_output=True, text=True, check=True).stdout
    h, v = [l.split('\t') for l in st_out.strip().split('\n')]
    d = dict(zip(h, v))

    fl = open(os.path.join(C.WORK, 'depth', s + '.flagstat.txt')).read()
    pct = float(next(l for l in fl.split('\n') if ' mapped (' in l).split('(')[1].split('%')[0])

    dep = [float(l.split()[1]) for l in open(os.path.join(C.WORK, 'depth', s + '.chrom1kb.txt'))]
    med = statistics.median(dep) if dep else 0.0

    rows.append({'sample': s, 'day': C.DAY.get(s, ''), 'rep': C.REP.get(s, ''),
                 'temp_C': C.TEMP.get(s, ''), 'condition': C.CONDITION.get(s, ''),
                 'reads': d['num_seqs'], 'bases': d['sum_len'], 'mean_len': d['avg_len'],
                 'N50': d['N50'], 'max_len': d['max_len'], 'Q20_pct': d.get('Q20(%)', ''),
                 'mean_qual': d.get('AvgQual', ''), 'GC_pct': d.get('GC(%)', ''),
                 'pct_mapped': pct, 'median_chrom_depth': round(med, 1)})
    print('%-8s reads=%-9s depth=%.1f mapped=%.2f%%' % (s, d['num_seqs'], med, pct))

out = os.path.join(C.ROOT, 'read_stats.csv')
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print('wrote', out)
