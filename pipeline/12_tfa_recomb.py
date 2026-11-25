"""Frequency of lambda tail-fibre (tfa) alleles carrying Rac-prophage donor sequence."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
WORK = C.WORK
import subprocess, re, json, math, collections
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
S = C.SAMPLES
W0, W1, F = 3083700, 3084600, 150

def wilson(x, n, z=1.96):
    if not n: return (None, None)
    d = 1 + z*z/n; c = (x/n + z*z/(2*n))/d
    h = z*math.sqrt((x/n)*(1-x/n)/n + z*z/(4*n*n))/d
    return max(0., c-h), min(1., c+h)

res = {}
for s in S:
    rec = wt = 0; pos_c = collections.Counter()
    p = subprocess.run(['samtools', 'view', '-F', '0x904', '-q', '20',
                        C.bam(s), 'contig_1:%d-%d' % (W0-2500, W1+2500)],
                       capture_output=True, text=True).stdout
    for line in p.split('\n'):
        if not line: continue
        f = line.split('\t', 10)
        pos, cig = int(f[3]), f[5]
        if cig == '*': continue
        r, ins, clip = pos, [], []
        for i, (n, o) in enumerate(CIG.findall(cig)):
            n = int(n)
            if o == 'I':
                if n >= 250 and W0 <= r <= W1: ins.append(r)
            elif o in 'M=X': r += n
            elif o in 'DN': r += n
            elif o in 'SH':
                if n >= 250: clip.append(pos if i == 0 else r)
        end = r
        clipin = [c for c in clip if W0 <= c <= W1]
        spans, enters = pos <= W0-F and end >= W1+F, (pos <= W0-F or end >= W1+F)
        if ins: rec += 1; pos_c[min(ins)//100*100] += 1
        elif clipin and enters: rec += 1; pos_c[('clip', min(clipin)//100*100)] += 1
        elif spans: wt += 1
    n = rec + wt
    res[s] = {'rec': rec, 'wt': wt, 'n': n, 'f': rec/n if n else None,
              'ci': wilson(rec, n), 'junctions': len(pos_c)}
json.dump(res, open(os.path.join(C.OUT, 'tfa.json'), 'w'), indent=1)
print('%-5s %5s %5s %5s  %-22s %s' % ('samp','rec','wt','n','f_recombinant [95% CI]','junctions'))
for s in S:
    v = res[s]
    if not v['n']: print('%-5s     - no spanning reads' % s); continue
    print('%-5s %5d %5d %5d  %.3f [%.3f-%.3f]      %d' %
          (s, v['rec'], v['wt'], v['n'], v['f'], *v['ci'], v['junctions']))
