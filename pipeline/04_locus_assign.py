#!/usr/bin/env python
"""Direct read assignment at the cassette locus.

Every read that covers enough of the variant-discriminating interior is
assigned to exactly ONE of the four alleles.  This is the method of record for
marker frequency; depth ratios at this locus are not trustworthy (coverage
shadow around the mask plus an induction-scaled gradient across the prophage).

The discriminating window deliberately stops short of the shared lysis block at
the 3' end of the interior, so a read is assigned on cassette sequence alone.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import subprocess, re, json, random, statistics as st

WORK = C.WORK
S = C.SAMPLES
V = C.VARIANTS
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
B = C.BOOT
SHARED_TAIL = 800          # bp of interior shared with the lysis block, excluded
random.seed(7)

REC = json.load(open(os.path.join(WORK, 'locus_ref.json')))
AL = C.ANCHOR_L[1] - C.ANCHOR_L[0] + 1
DISC = {v: (AL + 1, AL + REC[v]['interior_len'] - SHARED_TAIL) for v in V}

os.makedirs(os.path.join(WORK, 'locus'), exist_ok=True)
MMI = os.path.join(WORK, 'locus', 'ref.mmi')
subprocess.run(['minimap2', '-d', MMI, C.LOCUS_FASTA], capture_output=True, check=True)


def aligned_in(cig, pos, a, b):
    """Aligned (M/=/X) bases of this alignment falling inside [a, b]."""
    n = 0
    for num, op in CIG.findall(cig):
        num = int(num)
        if op in 'M=X':
            n += max(0, min(pos + num - 1, b) - max(pos, a) + 1)
            pos += num
        elif op in 'DN':
            pos += num
    return n


cnt, mapq = {}, {}
for s in S:
    fq = os.path.join(WORK, 'locus', s + '.fq')
    if not os.path.exists(fq):
        ns = os.path.join(WORK, 'locus', s + '.ns.bam')
        regions = '%s:%d-%d %s' % (C.CONTIG, C.CASSETTE_WINDOW[0], C.CASSETTE_WINDOW[1],
                                   ' '.join(C.ALLELES))
        subprocess.run(
            'samtools view -b -F 0x100 %s %s | samtools sort -n -o %s - '
            '&& samtools fastq -0 %s -s /dev/null %s 2>/dev/null'
            % (C._q(C.bam(s)), regions, C._q(ns), C._q(fq), C._q(ns)),
            shell=True, check=True)
        os.remove(ns)
    p = subprocess.run(
        'minimap2 -ax map-ont --secondary=no -t %d %s %s 2>/dev/null | samtools view -F 0x904 -'
        % (C.THREADS, C._q(MMI), C._q(fq)), shell=True, capture_output=True, text=True)
    c = {v: 0 for v in V}
    q = []
    for ln in p.stdout.splitlines():
        f = ln.split('\t', 6)
        v = f[2].replace('var_', '')
        if v not in c:
            continue
        a, b = DISC[v]
        if aligned_in(f[5], int(f[3]), a, b) >= C.MIN_DISC:
            c[v] += 1
            q.append(int(f[4]))
    cnt[s] = c
    mapq[s] = (st.median(q) if q else 0, sum(1 for x in q if x >= C.MAPQ) / max(1, len(q)))

# multinomial bootstrap on assigned reads
out = {}
for s in S:
    pool = [v for v in V for _ in range(cnt[s][v])]
    n = len(pool)
    reps = [random.choices(pool, k=n) for _ in range(B)] if n else []
    out[s] = {'n': n, 'mapq_med': mapq[s][0], 'mapq20': mapq[s][1]}
    for v in V:
        f = cnt[s][v] / n if n else 0.0
        bs = sorted(r.count(v) / n for r in reps) if n else [0.0]
        out[s][v] = {'n': cnt[s][v], 'f': f,
                     'lo': bs[int(.16 * len(bs))], 'hi': bs[int(.84 * len(bs))]}
json.dump(out, open(os.path.join(C.OUT, 'locus.json'), 'w'), indent=1)

# specificity check: a pure control must return its own allele and nothing else
for s, allele in C.PURE_CONTROLS.items():
    if s not in out:
        continue
    k = allele.replace('allele_', '')
    print('control %-6s %-4s f=%.4f  (off-target %.4f)'
          % (s, k, out[s][k]['f'], 1 - out[s][k]['f']))

print("sample  reads  mapQmed f>=20   " + "".join("%-18s" % ('f_' + v) for v in V))
for s in S:
    print("%-6s %6d  %4d  %5.2f   " % (s, out[s]['n'], out[s]['mapq_med'], out[s]['mapq20']) +
          "".join("%6.3f[%5.3f-%5.3f] " % (out[s][v]['f'], out[s][v]['lo'], out[s][v]['hi'])
                  for v in V))
