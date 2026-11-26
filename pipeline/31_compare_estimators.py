#!/usr/bin/env python
"""Side-by-side: stage 05 confusion-matrix unmixing vs stage 05b raw rates.

Neither is assumed correct.  This writes one row per library with both
estimates, their difference, and the interval widths, so the choice between
them is made on numbers rather than on argument.

Reads handoff/final.json (05) and handoff/raw_rates.json (05b).
Writes <OUT_TABLES>/estimator_comparison.csv and prints a summary.

What to look at:
  d_cmR_kanR     the only quantity that needs either model.  Renormalisation
                 cancels in a ratio, so a non-zero difference here comes from
                 the per-row matrix disagreeing with the pooled b, not from
                 bookkeeping.
  ci_ratio       width of the 05b interval over the 05 interval.  Expected > 1
                 wherever 05's fixed-M bootstrap was understating error.
  oob_05         stage 05 returned a fraction outside [0, 1].
  clamped_05     stage 05 returned exactly 0.0 where 05b returns > 0 -- the
                 max(0, .) clamp having hidden a negative point estimate.
"""
import os, sys, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

A = json.load(open(os.path.join(C.OUT, 'final.json')))          # stage 05
Braw = json.load(open(os.path.join(C.OUT, 'raw_rates.json')))    # stage 05b
Bs = Braw['samples']

TAB = getattr(C, 'OUT_TABLES', os.path.join(C.ROOT, 'results', 'tables'))
os.makedirs(TAB, exist_ok=True)
dest = os.path.join(TAB, 'estimator_comparison.csv')

F = ['p_Ia', 'p_IIa', 'p_cat']
rows = []
for s in C.SAMPLES:
    a, b = A[s], Bs[s]
    r = {'sample': s, 'rep': C.REP.get(s, ''), 'day': C.DAY.get(s, ''),
         'temp_C': C.TEMP.get(s, ''), 'condition': C.CONDITION.get(s, ''),
         'n_locus': b['n_locus'], 'b': round(Braw['b'], 4)}
    for k in F:
        r['%s_05' % k] = round(a[k]['f'], 4)
        r['%s_05b' % k] = round(b[k]['f'], 4)
        r['d_%s' % k] = round(b[k]['f'] - a[k]['f'], 4)
        wa = a[k]['hi'] - a[k]['lo']
        wb = b[k]['hi'] - b[k]['lo']
        r['ciw_%s_05' % k] = round(wa, 4)
        r['ciw_%s_05b' % k] = round(wb, 4)
        r['ci_ratio_%s' % k] = round(wb / wa, 3) if wa > 0 else ''
    ra, rb = a['cat_kan']['r'], b['cat_kan']['r']
    r['cmR_kanR_05'] = round(ra, 4) if ra not in (float('inf'),) else 'inf'
    r['cmR_kanR_05b'] = round(rb, 4) if rb not in (float('inf'),) else 'inf'
    r['d_cmR_kanR'] = (round(rb - ra, 4)
                       if float('inf') not in (ra, rb) else '')
    r['rel_d_cmR_kanR'] = (round((rb - ra) / ra, 4)
                           if float('inf') not in (ra, rb) and ra > 0 else '')
    r['cat_vs_Ia_05b'] = round(b['cat_vs_Ia']['f'], 4)
    r['cat_vs_Ia_n'] = b['cat_vs_Ia']['n']
    r['bor_excess_f'] = round(b['bor_excess_f']['f'], 4)
    r['bor_excess_lo'] = round(b['bor_excess_f']['lo'], 4)
    r['bor_excess_hi'] = round(b['bor_excess_f']['hi'], 4)
    r['sum_05'] = round(sum(a[k]['f'] for k in F), 4)
    r['sum_05b'] = round(sum(b[k]['f'] for k in F), 4)
    r['oob_05'] = ';'.join(k for k in F if a[k]['f'] < 0 or a[k]['f'] > 1) or ''
    r['clamped_05'] = ';'.join(k for k in F
                               if a[k]['f'] == 0.0 and b[k]['f'] > 0.002) or ''
    rows.append(r)

with open(dest, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print('wrote %s  (%d libraries)' % (dest, len(rows)))

# ---- summary -------------------------------------------------------------
def col(k):
    return [r[k] for r in rows if isinstance(r.get(k), float)]

print('\nb = %.4f  [%.4f-%.4f]  from %s'
      % (Braw['b'], Braw['b_lo'], Braw['b_hi'], '+'.join(Braw['b_controls'])))

d = [abs(x) for x in col('d_cmR_kanR')]
if d:
    print('\n|d cmR:kanR|  max=%.4f  median=%.4f' % (max(d), sorted(d)[len(d) // 2]))
rel = [abs(x) for x in col('rel_d_cmR_kanR')]
if rel:
    print('relative      max=%.1f%%  median=%.1f%%'
          % (100 * max(rel), 100 * sorted(rel)[len(rel) // 2]))

print('\nlargest cmR:kanR movements:')
ranked = sorted((r for r in rows if isinstance(r.get('rel_d_cmR_kanR'), float)),
                key=lambda r: -abs(r['rel_d_cmR_kanR']))
for r in ranked[:8]:
    print('  %-7s %-26s 05=%-8s 05b=%-8s  %+.1f%%'
          % (r['sample'], str(r['condition'])[:26], r['cmR_kanR_05'],
             r['cmR_kanR_05b'], 100 * r['rel_d_cmR_kanR']))

for k in F:
    cr = col('ci_ratio_%s' % k)
    if cr:
        print('\nCI width 05b/05 for %-5s median=%.2f  max=%.2f  (>1 = 05 was '
              'understating)' % (k, sorted(cr)[len(cr) // 2], max(cr)))

oob = [r for r in rows if r['oob_05']]
if oob:
    print('\nstage 05 out of [0,1]:')
    for r in oob:
        print('  %-7s %s  p_IIa_05=%s -> 05b=%s'
              % (r['sample'], r['oob_05'], r['p_IIa_05'], r['p_IIa_05b']))

cl = [r for r in rows if r['clamped_05']]
if cl:
    print('\nstage 05 clamped to exactly 0 where 05b is positive:')
    for r in cl:
        print('  %-7s %s' % (r['sample'], r['clamped_05']))

ex = sorted((r for r in rows if isinstance(r.get('bor_excess_f'), float)),
            key=lambda r: -r['bor_excess_f'])
print('\nbor in excess of own Delta-lysis population (top 5; '
      'interval excluding 0 = foreign bor):')
for r in ex[:5]:
    flag = '  <-- excludes 0' if r['bor_excess_lo'] > 0 else ''
    print('  %-7s %+.4f [%+.4f-%+.4f] n=%-5s %s%s'
          % (r['sample'], r['bor_excess_f'], r['bor_excess_lo'],
             r['bor_excess_hi'], r['n_locus'], str(r['condition'])[:24], flag))
