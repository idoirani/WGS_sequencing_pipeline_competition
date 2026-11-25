"""Per-junction IS family calls at lamB and malT, from insert length and donor locus."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
WORK = C.WORK
import subprocess, re, json, collections
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
S = C.SAMPLES
LOCI = {'lamB': C.GENES['lamB']['span'], 'malT': C.GENES['malT']['span']}
# IS family by element length (bp), from blast-confirmed chromosomal donors
FAM = [('IS1', 768, 40, 'IS1 (768 bp; donors 299.5/1904.5/3912.0 kb)'),
       ('IS5', 1195, 60, 'IS5 (1195 bp; donors 230.5/1780.8 kb)'),
       ('IS2', 1331, 70, 'IS2 (1331 bp; donor 4074.7 kb)')]
def fam(L):
    for n, c, tol, _ in FAM:
        if abs(L - c) <= tol: return n
    return 'other'

rows = []
for s in S:
    for g, (lo, hi) in LOCI.items():
        p = subprocess.run(['samtools', 'view', '-F', '0x904', '-q', '20',
                            C.bam(s), 'contig_1:%d-%d' % (lo, hi)],
                           capture_output=True, text=True).stdout
        ev = collections.defaultdict(list)
        for line in p.split('\n'):
            if not line: continue
            f = line.split('\t', 10)
            pos, cig = int(f[3]), f[5]
            if cig == '*': continue
            r = pos
            for n, o in CIG.findall(cig):
                n = int(n)
                if o == 'I':
                    if n >= 200 and lo <= r <= hi: ev[r // 80 * 80].append(n)
                elif o in 'M=X': r += n
                elif o in 'DN': r += n
        for b, ls in ev.items():
            if len(ls) < 2: continue
            med = sorted(ls)[len(ls) // 2]
            rows.append({'sample': s, 'gene': g, 'pos': b + 40, 'reads': len(ls),
                         'len': med, 'family': fam(med)})
rows.sort(key=lambda d: (d['gene'], d['sample'], d['pos']))
json.dump(rows, open(os.path.join(C.OUT, 'junctions.json'), 'w'), indent=1)
byfam = collections.Counter((r['gene'], r['family']) for r in rows)
print('distinct junctions (>=2 reads):', len(rows))
for k, v in sorted(byfam.items()): print('  %-6s %-6s %d' % (k[0], k[1], v))
print('\n%-5s %-5s %-9s %5s %6s %s' % ('samp','gene','pos','reads','len','family'))
for r in rows: print('%-5s %-5s %-9d %5d %6d %s' % (r['sample'], r['gene'], r['pos'], r['reads'], r['len'], r['family']))
