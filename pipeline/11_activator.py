"""Disrupted-receptor fraction: each read gets ONE class, one shared denominator per locus."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
WORK = C.WORK
import subprocess, re, json, collections
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
S = C.SAMPLES
WIN = {'lamB': (4325100, 4326040), 'malT': (330450, 330800)}
F = 150                                  # required flank outside the window

res = {}
for s in S:
    for g, (w0, w1) in WIN.items():
        p = subprocess.run(['samtools', 'view', '-F', '0x904', '-q', '20',
                            C.bam(s),
                            'contig_1:%d-%d' % (w0 - 2000, w1 + 2000)],
                           capture_output=True, text=True).stdout
        nref = nins = 0
        alleles = collections.Counter()
        for line in p.split('\n'):
            if not line: continue
            f = line.split('\t', 10)
            pos, cig = int(f[3]), f[5]
            if cig == '*': continue
            ref, ins, clip = pos, [], []
            ops = CIG.findall(cig)
            for i, (n, o) in enumerate(ops):
                n = int(n)
                if o == 'I':
                    if n >= 200: ins.append((ref, n))
                elif o in 'M=X': ref += n
                elif o in 'DN': ref += n
                elif o in 'SH':
                    if n >= 200: clip.append(pos if i == 0 else ref)
            end = ref
            inwin  = [(ip, il) for ip, il in ins if w0 <= ip <= w1]
            clipin = [cp for cp in clip if w0 <= cp <= w1]
            spans  = pos <= w0 - F and end >= w1 + F
            enters = (pos <= w0 - F) or (end >= w1 + F)
            if inwin:                                   # carries the element outright
                nins += 1; alleles[min(inwin)[0] // 60 * 60] += 1
            elif clipin and enters:                     # diverges inside the window
                nins += 1; alleles[('clip', min(clipin) // 60 * 60)] += 1
            elif spans:
                nref += 1
        tot = nref + nins
        res['%s|%s' % (s, g)] = {'ins': nins, 'ref': nref, 'n': tot,
                                 'f': nins / tot if tot else None,
                                 'lo': None, 'hi': None,
                                 'alleles': {str(k): v for k, v in alleles.items()}}
# Wilson 95% CI
import math
for k, v in res.items():
    n, x = v['n'], v['ins']
    if n:
        z = 1.96; ph = x / n; d = 1 + z * z / n
        c = (ph + z * z / (2 * n)) / d
        h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
        v['lo'], v['hi'] = max(0, c - h), min(1, c + h)
json.dump(res, open(os.path.join(C.OUT, 'recfreq2.json'), 'w'), indent=1)
print('%-5s %-5s %6s %6s %6s  %-22s %s' % ('samp','gene','ins','ref','n','f [95% CI]','distinct junctions'))
for g in WIN:
    for s in S:
        v = res['%s|%s' % (s, g)]
        if not v['n']: continue
        print('%-5s %-5s %6d %6d %6d  %.3f [%.3f-%.3f]   %d' %
              (s, g, v['ins'], v['ref'], v['n'], v['f'], v['lo'], v['hi'], len(v['alleles'])))
