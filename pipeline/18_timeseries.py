#!/usr/bin/env python
"""Per-condition time series: one row per condition group per passage day.

Metrics, all measured from reads:
  prophage_copies   median prophage depth / genome background.  At 37 C after
                    induction this counts replicating phage DNA as well as
                    integrated copies, so it is phage genome equivalents per
                    chromosome, not integrated prophage number.
  cmR_over_kanR     p_cat / (p_Ia + p_IIa), from unmixed marker fractions
  f_cmR             the same quantity as a fraction -- a ratio is unstable when
                    its denominator approaches zero, so both are reported
  lysis_minus_over_plus, f_lysis_minus   from S-R-Rz junction reads
  lamB_f_disrupted, malT_f_disrupted     from 16_lesions.py

Replicates are averaged only across flasks that exist for that group and day;
n_flasks is on every row.  With n <= 2 nothing is tested, and no row carries a
standard error implying otherwise -- the spread across the flasks present is
reported as min/max instead.
"""
import os, sys, csv, json, statistics as st
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

F = json.load(open(os.path.join(C.OUT, 'final.json')))
B = json.load(open(os.path.join(C.OUT, 'boot.json')))
LES = json.load(open(os.path.join(C.OUT, 'lesions.json')))['rows']
LESD = {(r['sample'], r['gene']): r for r in LES}
CON = {}
p = os.path.join(C.ROOT, 'genotype_concordance.csv')
if os.path.exists(p):
    CON = {r['sample']: r for r in csv.DictReader(open(p))}


def ratio(num, den):
    return num / den if den > 0 else ''


per = []
for s in C.SAMPLES:
    f = F[s]
    kan = f['p_Ia']['f'] + f['p_IIa']['f']
    cat = f['p_cat']['f']
    lp = f['lysis']['f']
    lam = LESD.get((s, 'lamB'), {})
    mal = LESD.get((s, 'malT'), {})
    per.append({
        'sample': s, 'rep': C.REP.get(s, ''), 'day': C.DAY.get(s, ''),
        'temp_C': C.TEMP.get(s, ''), 'group': C.GROUP.get(s, ''),
        'condition': C.CONDITION.get(s, ''),
        'concordance': CON.get(s, {}).get('verdict', ''),
        'n_locus_reads': f['n_locus'], 'n_junction_reads': f['lysis']['n'],
        'prophage_copies': round(B['cn'][s]['cn'], 3),
        'prophage_null_lo': round(B['cn'][s]['null_lo'], 3),
        'prophage_null_hi': round(B['cn'][s]['null_hi'], 3),
        'p_Ia': round(f['p_Ia']['f'], 4), 'p_IIa': round(f['p_IIa']['f'], 4),
        'p_cat': round(cat, 4),
        'f_cmR': round(ratio(cat, cat + kan), 4) if (cat + kan) else '',
        'cmR_over_kanR': round(ratio(cat, kan), 4) if kan > 0 else '',
        'f_lysis_minus': round(1 - lp, 4) if f['lysis']['n'] else '',
        'lysis_minus_over_plus': round(ratio(1 - lp, lp), 4) if lp > 0 else '',
        'lamB_n_reads': lam.get('n_reads', ''), 'lamB_f_disrupted': lam.get('f_disrupted', ''),
        'malT_n_reads': mal.get('n_reads', ''), 'malT_f_disrupted': mal.get('f_disrupted', ''),
    })

out1 = os.path.join(C.ROOT, 'timeseries_per_sample.csv')
with open(out1, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(per[0])); w.writeheader()
    w.writerows(sorted(per, key=lambda r: (r['group'], r['day'], r['rep'])))

METRICS = ['prophage_copies', 'f_cmR', 'cmR_over_kanR', 'f_lysis_minus',
           'lysis_minus_over_plus', 'lamB_f_disrupted', 'malT_f_disrupted']

bygroup = {}
for r in per:
    bygroup.setdefault((r['group'], r['day']), []).append(r)

agg = []
for (g, d), rs in sorted(bygroup.items()):
    row = {'group': g, 'day': d, 'temp_C': rs[0]['temp_C'],
           'condition': rs[0]['condition'], 'n_flasks': len(rs),
           'flasks': '+'.join(r['sample'] for r in rs),
           'flagged': '+'.join(r['sample'] for r in rs if r['concordance'] == 'FLAG')}
    for m in METRICS:
        v = [float(r[m]) for r in rs if r[m] not in ('', None)]
        row[m] = round(st.mean(v), 4) if v else ''
        row[m + '_min'] = round(min(v), 4) if v else ''
        row[m + '_max'] = round(max(v), 4) if v else ''
        row[m + '_n'] = len(v)
    agg.append(row)

out2 = os.path.join(C.ROOT, 'timeseries_by_condition.csv')
with open(out2, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(agg[0])); w.writeheader(); w.writerows(agg)

groups = sorted({r['group'] for r in per})
print('%d condition groups x up to 4 days -> %d aggregated rows' % (len(groups), len(agg)))
print('\n%-34s %-4s %-3s %8s %7s %8s %8s %8s' %
      ('group', 'day', 'n', 'copies', 'f_cmR', 'f_lys-', 'lamB', 'malT'))
for r in agg:
    print('%-34s %-4s %-3d %8s %7s %8s %8s %8s' %
          (r['group'][:34], r['day'], r['n_flasks'], r['prophage_copies'], r['f_cmR'],
           r['f_lysis_minus'], r['lamB_f_disrupted'], r['malT_f_disrupted']))
print('\nwrote %s, %s' % (out1, out2))
