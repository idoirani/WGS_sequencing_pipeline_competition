"""Genome-wide structural-variant scan from CIGARs: large I/D plus soft-clip clusters."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
WORK = C.WORK
import subprocess, re, json, collections

S = C.SAMPLES
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
MINSV, MINCLIP, WIN, MINREADS = 50, 200, 60, 4
MASK = (3100000, 3115000)            # engineered cassette locus - excluded

ev = collections.defaultdict(list)   # (kind,bin) -> [(sample, pos, len)]
tot = collections.Counter()
for s in S:
    p = subprocess.run(['samtools', 'view', '-F', '0x904', '-q', '20',
                        C.bam(s), C.CONTIG],
                       capture_output=True, text=True).stdout
    for line in p.split('\n'):
        if not line: continue
        f = line.split('\t', 6)
        pos, cig = int(f[3]), f[5]
        if cig == '*': continue
        tot[s] += 1
        ops = CIG.findall(cig)
        ref = pos
        for i, (n, o) in enumerate(ops):
            n = int(n)
            if o in 'MDN=X':
                if o in 'DN' and n >= MINSV and not (MASK[0] < ref < MASK[1]):
                    ev[('DEL', ref // WIN)].append((s, ref, n))
                ref += n
            elif o == 'I' and n >= MINSV and not (MASK[0] < ref < MASK[1]):
                ev[('INS', ref // WIN)].append((s, ref, n))
            elif o in 'SH' and n >= MINCLIP:
                b = pos if i == 0 else ref
                if not (MASK[0] < b < MASK[1]):
                    ev[('CLIP', b // WIN)].append((s, b, n))

out = []
for (kind, b), lst in ev.items():
    if len(lst) < MINREADS: continue
    per = collections.Counter(x[0] for x in lst)
    out.append({'kind': kind, 'pos': int(sum(x[1] for x in lst) / len(lst)),
                'n': len(lst), 'med_len': sorted(x[2] for x in lst)[len(lst) // 2],
                'samples': dict(per)})
out.sort(key=lambda d: -d['n'])
json.dump({'events': out, 'reads': dict(tot)}, open(os.path.join(C.OUT, 'sv.json'), 'w'), indent=1)
print('reads', sum(tot.values()), '| clusters >=%d reads:' % MINREADS, len(out))
for d in out[:28]:
    print('%-5s %9d  n=%-4d len=%-6d %s' % (d['kind'], d['pos'], d['n'], d['med_len'],
          ' '.join('%s:%d' % kv for kv in sorted(d['samples'].items(),
                                                 key=lambda k: -k[1])[:6])))
