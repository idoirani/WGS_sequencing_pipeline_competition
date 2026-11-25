#!/usr/bin/env python
"""Does each library's genotype match the condition it is labelled with?

Two independent observations per library -- which marker cassettes are present,
and the lysis genotype read at a different locus -- are compared with what the
lab's own condition label says should be in the flask.  A library that
disagrees is FLAGGED and reported.  It is never silently relabelled: a
disagreement can be a plating error, a pipetting error, or a transposed
barcode, and only the lab can say which.

CONDITION_MARKERS in config_inputs.py is the lab fact this stage tests against
(which strains go into which mixture); it is written out explicitly so the test
is not circular with the data.
"""
import os, sys, re, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import config_inputs as I

PRESENT = 0.05        # a marker counts as present above this unmixed fraction
LYSIS_TOL = 0.05      # tolerance on a lysis fraction that should be 0 or 1

F = json.load(open(os.path.join(C.OUT, 'final.json')))
KEY = {'Ia': 'p_Ia', 'IIa': 'p_IIa', 'cat': 'p_cat'}


def expected(cond):
    for pat, markers in I.CONDITION_MARKERS:
        if re.search(pat, cond, re.I):
            return set(markers)
    return None


rows = []
for s in C.SAMPLES:
    cond = C.CONDITION.get(s, '')
    exp = expected(cond)
    if exp is None:
        print('  ! no marker expectation matches condition %r (%s)' % (cond, s))
        continue
    obs = {m for m in KEY if F[s][KEY[m]]['f'] >= PRESENT}
    missing, extra = exp - obs, obs - exp
    # lysis expectation: pure lysis- is 0, anything with no IIa is 1, mixtures free
    lyf = F[s]['lysis']['f']
    if exp == {'IIa'}:
        lys_ok = lyf <= LYSIS_TOL
        lys_exp = '0'
    elif 'IIa' not in exp:
        lys_ok = lyf >= 1 - LYSIS_TOL
        lys_exp = '1'
    else:
        lys_ok, lys_exp = True, 'mixed'
    verdict = 'OK' if (not missing and not extra and lys_ok) else 'FLAG'
    rows.append({'sample': s, 'rep': C.REP.get(s, ''), 'day': C.DAY.get(s, ''),
                 'temp_C': C.TEMP.get(s, ''), 'condition': cond,
                 'expected_markers': '+'.join(sorted(exp)),
                 'observed_markers': '+'.join(sorted(obs)) or '-',
                 'p_Ia': round(F[s]['p_Ia']['f'], 4),
                 'p_IIa': round(F[s]['p_IIa']['f'], 4),
                 'p_cat': round(F[s]['p_cat']['f'], 4),
                 'n_locus_reads': F[s]['n_locus'],
                 'f_lysis_plus': round(lyf, 4), 'n_junction_reads': F[s]['lysis']['n'],
                 'lysis_expected': lys_exp,
                 'missing_markers': '+'.join(sorted(missing)) or '',
                 'unexpected_markers': '+'.join(sorted(extra)) or '',
                 'lysis_concordant': 'yes' if lys_ok else 'no',
                 'verdict': verdict})

out = os.path.join(C.ROOT, 'genotype_concordance.csv')
with open(out, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

flag = [r for r in rows if r['verdict'] == 'FLAG']
print('concordance: %d of %d libraries agree with their condition label' %
      (len(rows) - len(flag), len(rows)))
if flag:
    print('\nFLAGGED (reported, not relabelled):')
    print('  %-7s %-3s %-3s %-28s %-9s %-9s %s' %
          ('sample', 'rep', 'day', 'condition', 'expected', 'observed', 'why'))
    for r in flag:
        why = []
        if r['missing_markers']:
            why.append('missing ' + r['missing_markers'])
        if r['unexpected_markers']:
            why.append('unexpected ' + r['unexpected_markers'])
        if r['lysis_concordant'] == 'no':
            why.append('lysis+ %.3f, expected %s' % (r['f_lysis_plus'], r['lysis_expected']))
        print('  %-7s %-3s %-3s %-28s %-9s %-9s %s' %
              (r['sample'], r['rep'], r['day'], r['condition'][:28],
               r['expected_markers'], r['observed_markers'], '; '.join(why)))
print('\nwrote', out)
