"""One read -> one lamB allele class, single denominator spanning all three lesion positions."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
WORK = C.WORK
import subprocess, re, json, math, collections
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
S = C.SAMPLES
W0, W1, F = 4325100, 4326040, 150
FS, NS = 4325749, 4325797            # frameshift homopolymer / nonsense SNV

def wilson(x, n, z=1.96):
    if not n: return (None, None)
    d = 1 + z * z / n; c = (x / n + z * z / (2 * n)) / d
    h = z * math.sqrt((x / n) * (1 - x / n) / n + z * z / (4 * n * n)) / d
    return max(0., c - h), min(1., c + h)

res = {}
for s in S:
    cnt = collections.Counter(); junc = collections.Counter()
    p = subprocess.run(['samtools', 'view', '-F', '0x904', '-q', '20',
                        C.bam(s), 'contig_1:%d-%d' % (W0 - 2500, W1 + 2500)],
                       capture_output=True, text=True).stdout
    for line in p.split('\n'):
        if not line: continue
        f = line.split('\t', 10)
        pos, cig, seq = int(f[3]), f[5], f[9]
        if cig == '*': continue
        r, q = pos, 0
        isins, fs, clip, cov_fs, cov_ns, ns = [], False, [], False, False, False
        for i, (n, o) in enumerate(CIG.findall(cig)):
            n = int(n)
            if o == 'I':
                if n >= 200 and W0 <= r <= W1: isins.append(r)
                if 1 <= n <= 5 and abs(r - FS) <= 6: fs = True
                q += n
            elif o in 'M=X':
                if r <= FS <= r + n:
                    cov_fs = True
                if r <= NS <= r + n:
                    cov_ns = True
                    if seq[q + NS - r].upper() == 'T': ns = True
                r += n; q += n
            elif o in 'DN':
                if n <= 5 and abs(r - FS) <= 6: fs = True
                r += n
            elif o in 'SH':
                if n >= 200: clip.append(pos if i == 0 else r)
                q += n if o == 'S' else 0
        end = r
        spans  = pos <= W0 - F and end >= W1 + F
        enters = pos <= W0 - F or end >= W1 + F
        clipin = [c for c in clip if W0 <= c <= W1]
        if isins:
            cnt['IS'] += 1; junc[min(isins) // 60 * 60] += 1
        elif clipin and enters:
            cnt['IS'] += 1; junc[('clip', min(clipin) // 60 * 60)] += 1
        elif fs and cov_fs:
            cnt['FS'] += 1
        elif ns and cov_ns:
            cnt['NS'] += 1
        elif spans:
            cnt['WT'] += 1
    n = sum(cnt.values()); dis = n - cnt['WT']
    res[s] = {'n': n, **{k: cnt[k] for k in ('IS', 'FS', 'NS', 'WT')},
              'f_dis': dis / n if n else None, 'ci': wilson(dis, n),
              'n_junctions': len(junc)}
json.dump(res, open(os.path.join(C.OUT, 'lamb_final.json'), 'w'), indent=1)
print('%-5s %5s %5s %5s %5s %5s  %-22s %s' % ('samp','n','IS','FS','NS','WT','f_disrupted [95% CI]','junctions'))
for s in S:
    v = res[s]
    print('%-5s %5d %5d %5d %5d %5d  %.3f [%.3f-%.3f]      %d' %
          (s, v['n'], v['IS'], v['FS'], v['NS'], v['WT'], v['f_dis'], *v['ci'], v['n_junctions']))
