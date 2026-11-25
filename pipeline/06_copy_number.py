#!/usr/bin/env python
"""Prophage copy number per chromosome, with a faux-prophage null.

Copy number is the median 1 kb depth inside the prophage divided by the
genome-wide background median.  The error bar is NOT a bootstrap of that ratio:
it is the spread of the SAME statistic computed on same-length windows of
ordinary chromosome ("faux prophages"), which is what tells you how much of an
apparent gain is just regional depth variation.

Also bootstraps the depth-based allele and lysis fractions, calibrated on the
pure controls.  Those depth fractions are diagnostic only -- stage 04/05 read
assignment is the method of record (see TOOLS.md).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import subprocess, re, json, random, statistics as st

WORK = C.WORK
S = C.SAMPLES
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
B = C.BOOT
random.seed(0)

META = json.load(open(os.path.join(WORK, 'ref_probe.json')))
REG = {k: (k, META[k]['cds'][0], META[k]['cds'][1]) for k in C.ALLELES}
REG['lysis'] = (C.CONTIG, C.REGION['lysis'][0], C.REGION['lysis'][1])
REG['denom'] = (C.CONTIG, C.REGION['denom'][0], C.REGION['denom'][1])


def read_overlaps(bam, contig, a, b):
    """Aligned (M/=/X) bases per read falling inside [a,b]. One entry per read."""
    out = []
    p = subprocess.run(['samtools', 'view', '-F', '0x904', bam, '%s:%d-%d' % (contig, a, b)],
                       capture_output=True, text=True, check=True)
    for ln in p.stdout.splitlines():
        f = ln.split('\t', 6)
        pos, n = int(f[3]), 0
        for num, op in CIG.findall(f[5]):
            num = int(num)
            if op in 'M=X':
                n += max(0, min(pos + num - 1, b) - max(pos, a) + 1)
                pos += num
            elif op in 'DN':
                pos += num
        if n:
            out.append(n)
    return out


def boot_depth(ov, L, B):
    if not ov:
        return [0.0] * B
    k = len(ov)
    return [sum(random.choices(ov, k=k)) / L for _ in range(B)]


# ---- 1. read-level bootstrap of every region in every sample -------------
obs, reps = {}, {}
for s in S:
    for r, (c, a, b) in REG.items():
        ov = read_overlaps(C.bam(s), c, a, b)
        L = b - a + 1
        obs[(s, r)] = sum(ov) / L
        reps[(s, r)] = boot_depth(ov, L, B)

# ---- 2. fractions, calibrated on the pure controls ----------------------
CTRL = {a: list(ss) for a, ss in C.CONTROLS_BY_ALLELE.items() if a in REG}
lysis_ctrl = sorted({s for a in ('allele_Ia', 'allele_cat')
                     for s in C.CONTROLS_BY_ALLELE.get(a, [])})
if lysis_ctrl:
    CTRL['lysis'] = lysis_ctrl


def ratio_reps(s, r):
    return [n / d if d > 0 else float('nan')
            for n, d in zip(reps[(s, r)], reps[(s, 'denom')])]


frac = {}
for r, ctrls in CTRL.items():
    ctrls = [c for c in ctrls if c in S]
    if not ctrls:
        continue
    cal_obs = st.mean(obs[(c, r)] / obs[(c, 'denom')] for c in ctrls)
    cal_reps = [st.mean(v) for v in zip(*[ratio_reps(c, r) for c in ctrls])]
    for s in S:
        o = (obs[(s, r)] / obs[(s, 'denom')]) / cal_obs if cal_obs else float('nan')
        rr = sorted(x / y for x, y in zip(ratio_reps(s, r), cal_reps))
        frac[(s, r)] = {'f': o, 'lo': rr[int(.16 * B)], 'hi': rr[int(.84 * B)],
                        'sd': st.pstdev(rr), 'depth': obs[(s, r)]}

# ---- 3. prophage copy number + faux-prophage null ------------------------
WL = (C.PROPHAGE[1] - C.PROPHAGE[0]) // 1000          # window length in 1 kb bins


def excluded(p):
    return (C.PROPHAGE[0] - 5000 <= p <= C.PROPHAGE[1] + 5000 or
            any(a <= p <= b for a, b in C.CRYPTIC))


cn = {}
for s in S:
    d = {}
    for ln in open(os.path.join(WORK, 'depth', s + '.chrom1kb.txt')):
        p, v = ln.split()
        d[int(p)] = float(v)
    pos = sorted(d)
    bg = st.median(d[p] for p in pos if not excluded(p) and d[p] > 0)
    pro = st.median(d[p] for p in pos if C.PROPHAGE[0] <= p <= C.PROPHAGE[1])
    # null: same-length windows of contiguous non-excluded bins elsewhere
    runs, cur = [], []
    for p in pos:
        if excluded(p) or d[p] <= 0:
            if len(cur) > WL:
                runs.append(cur)
            cur = []
        else:
            cur.append(p)
    if len(cur) > WL:
        runs.append(cur)
    faux = []
    for _ in range(B):
        r = random.choices(runs, weights=[len(x) - WL for x in runs])[0]
        i = random.randrange(len(r) - WL)
        faux.append(st.median(d[p] for p in r[i:i + WL]) / bg)
    faux.sort()
    cn[s] = {'cn': pro / bg, 'bg_depth': bg, 'null_mean': st.mean(faux),
             'null_sd': st.pstdev(faux), 'null_lo': faux[int(.16 * B)],
             'null_hi': faux[int(.84 * B)], 'n_windows': len(runs)}

json.dump({'frac': {'%s|%s' % (s, r): v for (s, r), v in frac.items()}, 'cn': cn,
           'obs_depth': {'%s|%s' % (s, r): v for (s, r), v in obs.items()}},
          open(os.path.join(C.OUT, 'boot.json'), 'w'), indent=1)
print('regions', len(REG), 'samples', len(S), 'B', B)
for r, ctrls in CTRL.items():
    ctrls = [c for c in ctrls if c in S]
    if ctrls:
        print('  calibrator %-12s %s' % (r, ' '.join(
            '%s %.4f' % (c, obs[(c, r)] / obs[(c, 'denom')]) for c in ctrls)))
print('%-6s %7s %9s %s' % ('sample', 'cn', 'bg_depth', 'null [16-84%]'))
for s in S:
    v = cn[s]
    print('%-6s %7.2f %9.1f  %.2f-%.2f' % (s, v['cn'], v['bg_depth'], v['null_lo'], v['null_hi']))
