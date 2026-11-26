#!/usr/bin/env python
"""Alternative marker estimator: raw read rates, one structural parameter.

Stage 05 inverts a 3x3 confusion matrix fitted per plate row.  This stage is a
replacement for the marker half of it (Part 2 only; the lysis junction in 05
Part 1 is independent of all of this and is not touched).  It exists because
the thing 05 calls an "IIa -> wild-type leak" is not a leak:

  Ia  and cat  REPLACE the bor-p78-p79 segment.
  IIa REPLACES S-R-Rz and leaves bor INTACT.

So a lysis-deficient genome contributes locus reads to TWO variants -- its own
IIa cassette and its own untouched bor -- while an Ia or cat genome contributes
only to its cassette.  Nothing is mis-assigned and nothing crossed between
wells.  The measured day-8 numbers say exactly this: bor_WT is 0.000-0.015 in
every library with no IIa, and 0.063-0.184 in every library with IIa.

One parameter follows:

  b = P(a Delta-lysis genome's locus read lands on intact bor, not on IIa)

which is read-length / window geometry, not a property of a plate row.  It is
estimated ONCE, pooled over the read counts of every pure IIa control, so it is
weighted by depth instead of averaging six fractions measured from 28-338
reads.  Genome-equivalent weights are then

  w_Ia  = n_Ia          w_cat = n_cat          w_IIa = n_IIa / (1 - b)

and the reported fractions are w / sum(w).  Consequences, versus stage 05:

  * every fraction is in [0, 1] by construction -- no max(0, .) clamp, so a
    true zero is estimated as zero instead of one-sidedly upward, and p_IIa
    cannot come back above 1 (stage 05 returns 1.0192 for v2D10);
  * b's own sampling error is resampled inside the bootstrap, so the intervals
    include control uncertainty.  Stage 05 holds M fixed across all B
    replicates and therefore understates its error bars;
  * cat:Ia is reported separately as a plain binomial.  Both cassettes sit at
    the SAME position behind byte-identical anchors, so that comparison is
    head-to-head and needs no model at all.

Diagnostic.  With b fixed, a library's expected bor count is b/(1-b) * n_IIa.
Any EXCESS over that is bor sequence not explained by its own Delta-lysis
population -- an unmarked/WT genome, or carryover.  Reported per library as
bor_excess_f with a bootstrap interval; it is the contamination channel that a
confusion matrix cannot express, because the matrix absorbs bor into a diagonal.

Writes handoff/raw_rates.json.  Stage 05 is left intact: run both, then
31_compare_estimators.py for the side-by-side.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import json, random, math

random.seed(13)
B = C.BOOT
CASS = ('Ia', 'IIa', 'cat')
L = json.load(open(os.path.join(C.OUT, 'locus.json')))


def wilson(k, n, z=1.0):
    """Wilson interval at ~68% (z=1) to match the bootstrap quantiles used here."""
    if n == 0:
        return (float('nan'), float('nan'), float('nan'))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (p, max(0.0, c - h), min(1.0, c + h))


# ---- b, pooled over the read counts of every pure IIa control -------------
# Only controls that are actually pure contribute: a control carrying another
# cassette above CTRL_MAX_OFF_CASSETTE is excluded here for the same reason
# stage 05 excludes it from its matrix.
iia_ctrl, n_iia, n_bor = [], 0, 0
for s in C.CONTROLS_BY_ALLELE.get('allele_IIa', []):
    if s not in L:
        continue
    off = sum(L[s][v]['n'] for v in CASS if v != 'IIa')
    tot = max(1, L[s]['n'])
    if off / tot > C.CTRL_MAX_OFF_CASSETTE:
        print('  ! pure IIa control %s excluded from b (%.1f%% other cassette)'
              % (s, 100 * off / tot))
        continue
    iia_ctrl.append(s)
    n_iia += L[s]['IIa']['n']
    n_bor += L[s]['bor']['n']

if not iia_ctrl:
    raise SystemExit('no usable pure IIa control -- b cannot be estimated. Add one '
                     '(control_allele=allele_IIa) or run stage 05 instead.')

b_hat = n_bor / (n_iia + n_bor)
# control read pool for resampling b
ctrl_pool = ['bor'] * n_bor + ['IIa'] * n_iia
print('b (intact-bor share of a Delta-lysis genome\'s locus reads) = %.4f'
      '   from %d reads over %s' % (b_hat, len(ctrl_pool), '+'.join(iia_ctrl)))
_, b_lo, b_hi = wilson(n_bor, n_iia + n_bor)
print('    68%% interval %.4f-%.4f' % (b_lo, b_hi))

# per-control value, printed so a genuine row effect would still be visible
print('    per control: ' + '  '.join(
    '%s=%.4f(n=%d)' % (s, L[s]['bor']['n'] / max(1, L[s]['bor']['n'] + L[s]['IIa']['n']),
                       L[s]['bor']['n'] + L[s]['IIa']['n']) for s in iia_ctrl))


def weights(n_Ia, n_IIa, n_cat, b):
    w = (n_Ia, n_IIa / (1 - b) if b < 1 else 0.0, n_cat)
    t = sum(w)
    return (w[0] / t, w[1] / t, w[2] / t) if t > 0 else (0.0, 0.0, 0.0)


out = {}
for s in C.SAMPLES:
    r = {v: L[s][v]['n'] for v in ('Ia', 'IIa', 'cat', 'bor')}
    n_tot = sum(r.values())
    pt = weights(r['Ia'], r['IIa'], r['cat'], b_hat)

    # two-level bootstrap: the library's reads AND the control reads for b,
    # so the intervals carry b's sampling error as well as the library's.
    pool = [v for v in ('Ia', 'IIa', 'cat', 'bor') for _ in range(r[v])]
    reps, exc = [], []
    for _ in range(B):
        bb = random.choices(ctrl_pool, k=len(ctrl_pool)).count('bor') / len(ctrl_pool)
        if n_tot:
            d = random.choices(pool, k=n_tot)
            reps.append(weights(d.count('Ia'), d.count('IIa'), d.count('cat'), bb))
            exc.append((d.count('bor') - bb / (1 - bb) * d.count('IIa')) / n_tot)
        else:
            reps.append((0.0, 0.0, 0.0))
            exc.append(0.0)
    exc.sort()

    e = {'n_locus': n_tot, 'b': b_hat, 'raw_counts': r}
    for i, v in enumerate(['p_Ia', 'p_IIa', 'p_cat']):
        q = sorted(x[i] for x in reps)
        e[v] = {'f': pt[i], 'lo': q[int(.16 * B)], 'hi': q[int(.84 * B)]}

    # cmR : kanR  -- kanR pools a bor-locus allele with an SRRz-locus allele,
    # so this is the one quantity that needs b at all.
    rr = sorted(x[2] / (x[0] + x[1]) for x in reps if (x[0] + x[1]) > 0)
    kan = pt[0] + pt[1]
    e['cat_kan'] = {'r': pt[2] / kan if kan > 0 else float('inf'),
                    'lo': rr[int(.16 * len(rr))] if len(rr) > 20 else float('nan'),
                    'hi': rr[int(.84 * len(rr))] if len(rr) > 20 else float('nan')}

    # cat : Ia  -- same position, identical anchors, no model
    f, lo, hi = wilson(r['cat'], r['cat'] + r['Ia'])
    e['cat_vs_Ia'] = {'n': r['cat'] + r['Ia'], 'f': f, 'lo': lo, 'hi': hi}

    # bor in excess of this library's own Delta-lysis population
    exp_bor = b_hat / (1 - b_hat) * r['IIa'] if b_hat < 1 else 0.0
    e['bor_excess_f'] = {'f': (r['bor'] - exp_bor) / n_tot if n_tot else 0.0,
                         'lo': exc[int(.16 * B)], 'hi': exc[int(.84 * B)]}
    out[s] = e

json.dump({'b': b_hat, 'b_lo': b_lo, 'b_hi': b_hi, 'b_controls': iia_ctrl,
           'samples': out}, open(os.path.join(C.OUT, 'raw_rates.json'), 'w'), indent=1)

print('\n%-7s %6s  %-16s %-16s %-16s %6s  %-16s %-16s'
      % ('sample', 'n_loc', 'p_Ia', 'p_IIa', 'p_cat', 'sum', 'cat:Ia', 'bor_excess'))
for s in C.SAMPLES:
    e = out[s]
    tot = e['p_Ia']['f'] + e['p_IIa']['f'] + e['p_cat']['f']
    print('%-7s %6d  ' % (s, e['n_locus']) +
          ''.join('%6.3f[%5.3f-%5.3f] ' % (e[v]['f'], e[v]['lo'], e[v]['hi'])
                  for v in ('p_Ia', 'p_IIa', 'p_cat')) +
          '%6.3f  ' % tot +
          '%6.3f[%5.3f-%5.3f] ' % (e['cat_vs_Ia']['f'], e['cat_vs_Ia']['lo'],
                                   e['cat_vs_Ia']['hi']) +
          '%6.3f[%6.3f-%6.3f]' % (e['bor_excess_f']['f'], e['bor_excess_f']['lo'],
                                  e['bor_excess_f']['hi']))

print('\nControls must return their own allele at 1.000 and bor_excess at 0.000:')
for s, al in C.PURE_CONTROLS.items():
    if s in out:
        k = al.replace('allele_', '')
        print('  %-7s %-4s f=%.4f  bor_excess=%+.4f [%+.4f-%+.4f]'
              % (s, k, out[s]['p_' + k]['f'], out[s]['bor_excess_f']['f'],
                 out[s]['bor_excess_f']['lo'], out[s]['bor_excess_f']['hi']))
