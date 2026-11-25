#!/usr/bin/env python
"""Cassette-independent lysis status, plus unmixing of the marker assignment.

Part 1. Every read spanning the 3' boundary of the lysis block votes lysis+
(endolysin R + holin S present) or lysis- (the IIa cassette occupying the same
position).  This is a DIFFERENT locus from the marker assay, so agreement
between the two is evidence that no recombinant has uncoupled marker from
lysis genotype.

Part 2. Near-identical alleles cross-assign.  The pure single-strain controls
named in the sample sheet are the rows of a confusion matrix, which is inverted
analytically to recover true fractions.  Without controls this cannot be done
and the raw fractions are biased.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import subprocess, re, json, random, statistics as st
from Bio import SeqIO

WORK = C.WORK
S = C.SAMPLES
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
B = C.BOOT
MIN_DISC = C.MIN_DISC
random.seed(11)


def genome(strain):
    p = os.path.join(C.ANCESTOR_DIR, strain, strain + '_1_reference.fasta')
    return {r.id: str(r.seq) for r in SeqIO.parse(p, 'fasta')}


g_plus = genome(C.REF_STRAIN)
ms, mc, m0, m1 = C.LYSIS_MINUS
g_minus = genome(ms)
AR = g_plus[C.CONTIG][C.LYSIS_ANCHOR[0] - 1:C.LYSIS_ANCHOR[1]]
ctx = {'plus':  g_plus[C.CONTIG][C.LYSIS_CONTEXT[0] - 1:C.LYSIS_CONTEXT[1]],
       'minus': g_minus[mc][m0 - 1:m1]}
ref = os.path.join(WORK, 'lysjx_ref.fa')
open(ref, 'w').write(''.join('>lys_%s\n%s\n' % (k, v + AR) for k, v in ctx.items()))
DISC = (1, len(ctx['plus']))

os.makedirs(os.path.join(WORK, 'locus'), exist_ok=True)
MMI = os.path.join(WORK, 'locus', 'lysjx.mmi')
subprocess.run(['minimap2', '-d', MMI, ref], capture_output=True, check=True)

cnt = {}
for s in S:
    fq = os.path.join(WORK, 'locus', s + '.fq')
    p = subprocess.run(
        'minimap2 -ax map-ont --secondary=no -t %d %s %s 2>/dev/null | samtools view -F 0x904 -'
        % (C.THREADS, C._q(MMI), C._q(fq)), shell=True, capture_output=True, text=True)
    c = {'plus': 0, 'minus': 0}
    for ln in p.stdout.splitlines():
        f = ln.split('\t', 6)
        k = f[2].replace('lys_', '')
        if k not in c:
            continue
        a, b = DISC
        pos, n = int(f[3]), 0
        for num, op in CIG.findall(f[5]):
            num = int(num)
            if op in 'M=X':
                n += max(0, min(pos + num - 1, b) - max(pos, a) + 1)
                pos += num
            elif op in 'DN':
                pos += num
        if n >= MIN_DISC:
            c[k] += 1
    cnt[s] = c

# ---- confusion matrix, matched to each library's own plate row ------------
# Cross-assignment is NOT a constant of the assay: the IIa->wild-type leak
# measured here ranges 0.055-0.139 between plate rows.  Pooling all controls
# would push that spread into every mixture as a systematic error, so each
# library is unmixed with the controls from its OWN rep and day whenever those
# exist (the plate layout puts the three pure controls in columns 1-3 of every
# row).  Where a row's control is missing or failed QC, the pooled mean is used
# and the substitution is reported rather than hidden.
L = json.load(open(os.path.join(C.OUT, 'locus.json')))
CASS = ('Ia', 'IIa', 'cat')

ctrl_off = {}
for s, al in C.PURE_CONTROLS.items():
    if s in L:
        a = al.replace('allele_', '')
        ctrl_off[s] = (a, sum(L[s][v]['f'] for v in CASS if v != a))
failed = {s for s, (a, off) in ctrl_off.items() if off > C.CTRL_MAX_OFF_CASSETTE}


def _rows(k, key=None):
    out = []
    for s in C.CONTROLS_BY_ALLELE.get('allele_' + k, []):
        if s not in L or s in failed:
            continue
        if key is None or (C.REP.get(s), C.DAY.get(s)) == key:
            out.append(s)
    return out


_cache = {}


def matrix_for(sample):
    key = (C.REP.get(sample), C.DAY.get(sample))
    if key in _cache:
        return _cache[key]
    M, src = {}, {}
    for k in CASS:
        local = _rows(k, key)
        use = local or _rows(k)
        if not use:
            raise SystemExit('no usable pure control for allele_%s -- add one to the '
                             'sample sheet (or check which controls failed QC)' % k)
        M[k] = {v: st.mean(L[s][v]['f'] for s in use) for v in C.VARIANTS}
        src[k] = ('matched' if local else 'pooled', use)
    M['det'] = M['Ia']['Ia'] * M['cat']['cat'] - M['cat']['Ia'] * M['Ia']['cat']
    _cache[key] = (M, src)
    return M, src


def unmix(M, o_Ia, o_IIa, o_cat):
    p_IIa = o_IIa / M['IIa']['IIa']
    p_Ia = (M['cat']['cat'] * o_Ia - M['cat']['Ia'] * o_cat) / M['det']
    p_cat = (M['Ia']['Ia'] * o_cat - M['Ia']['cat'] * o_Ia) / M['det']
    return max(0.0, p_Ia), max(0.0, p_IIa), max(0.0, p_cat)


if failed:
    print('CONTROL QC FAILED (excluded from every confusion matrix):')
    for s in sorted(failed):
        a, off = ctrl_off[s]
        comp = ' '.join('%s=%.3f' % (v, L[s][v]['f']) for v in C.VARIANTS)
        print('  %-7s declared pure %-4s but %.1f%% of reads carry another cassette  [%s]'
              % (s, a, 100 * off, comp))


out = {}
for s in S:
    M, src = matrix_for(s)
    n = L[s]['n']
    pool = [v for v in C.VARIANTS for _ in range(L[s][v]['n'])]
    reps = []
    for _ in range(B):
        r = random.choices(pool, k=n) if n else []
        reps.append(unmix(M, r.count('Ia') / n, r.count('IIa') / n, r.count('cat') / n)
                    if n else (0, 0, 0))
    pt = unmix(M, L[s]['Ia']['f'], L[s]['IIa']['f'], L[s]['cat']['f'])
    e = {}
    for i, v in enumerate(['p_Ia', 'p_IIa', 'p_cat']):
        bs = sorted(x[i] for x in reps)
        e[v] = {'f': pt[i], 'lo': bs[int(.16 * B)], 'hi': bs[int(.84 * B)]}
    # cmR : kanR on the unmixed scale
    rr = sorted(x[2] / (x[0] + x[1]) for x in reps if (x[0] + x[1]) > 0)
    e['cat_kan'] = {'r': pt[2] / (pt[0] + pt[1]) if (pt[0] + pt[1]) > 0 else float('inf'),
                    'lo': rr[int(.16 * len(rr))] if len(rr) > 20 else float('nan'),
                    'hi': rr[int(.84 * len(rr))] if len(rr) > 20 else float('nan')}
    # lysis junction: binomial
    cp, cm = cnt[s]['plus'], cnt[s]['minus']
    nn = cp + cm
    bs = sorted(sum(random.random() < cp / nn for _ in range(nn)) / nn
                for _ in range(B)) if nn else [0.0]
    e['lysis'] = {'n': nn, 'n_plus': cp, 'f': cp / nn if nn else float('nan'),
                  'lo': bs[int(.16 * len(bs))], 'hi': bs[int(.84 * len(bs))]}
    e['n_locus'] = n
    out[s] = e
json.dump(out, open(os.path.join(C.OUT, 'final.json'), 'w'), indent=1)

print('\nconfusion matrix per plate row:')
for key in sorted(_cache, key=lambda k: (str(k[0]), k[1] if k[1] is not None else 0)):
    M, src = _cache[key]
    print('  rep %s day %-2s  ' % (key[0], key[1]) +
          '  '.join('%s:%s(%s)' % (k, src[k][0], '+'.join(src[k][1])) for k in CASS))
print("sample  n_loc  p_Ia            p_IIa           p_cat           sum    cat:kan            "
      "n_jx  f_lysis+        pred(Ia+cat)")
for s in S:
    e = out[s]
    tot = e['p_Ia']['f'] + e['p_IIa']['f'] + e['p_cat']['f']
    pred = e['p_Ia']['f'] + e['p_cat']['f']
    print("%-6s %5d  %s %s %s %5.3f  %-18s %4d  %s %6.3f  %+.3f" % (
        s, e['n_locus'],
        *("%6.3f[%5.3f-%5.3f]" % (e[v]['f'], e[v]['lo'], e[v]['hi'])
          for v in ('p_Ia', 'p_IIa', 'p_cat')),
        tot, "%6.3f[%5.3f-%5.3f]" % (e['cat_kan']['r'], e['cat_kan']['lo'], e['cat_kan']['hi']),
        e['lysis']['n'], "%6.3f[%5.3f-%5.3f]" % (e['lysis']['f'], e['lysis']['lo'], e['lysis']['hi']),
        pred, e['lysis']['f'] - pred))
