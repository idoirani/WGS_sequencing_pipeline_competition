"""Pull inserted sequences >=200 bp from reads overlapping lamB / malT, with exact junctions."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
WORK = C.WORK
import subprocess, re, collections
CIG = re.compile(r'(\d+)([MIDNSHP=X])')
LOCI = {'lamB': C.GENES['lamB']['span'], 'malT': C.GENES['malT']['span']}
S = C.SAMPLES
seqs, junc = [], collections.Counter()
for s in S:
    for nm, (lo, hi) in LOCI.items():
        p = subprocess.run(['samtools', 'view', '-F', '0x904', '-q', '20',
                            C.bam(s), 'contig_1:%d-%d' % (lo, hi)],
                           capture_output=True, text=True).stdout
        for line in p.split('\n'):
            if not line: continue
            f = line.split('\t')
            pos, cig, seq = int(f[3]), f[5], f[9]
            if cig == '*': continue
            ref, q = pos, 0
            for n, o in CIG.findall(cig):
                n = int(n)
                if o == 'I':
                    if n >= 200 and lo <= ref <= hi:
                        junc[(s, nm, ref, n // 50 * 50)] += 1
                        seqs.append(('%s_%s_%d_%d' % (s, nm, ref, n), seq[q:q + n]))
                    q += n
                elif o in 'M=X': ref += n; q += n
                elif o in 'DN': ref += n
                elif o == 'S': q += n
with open(WORK + '/ins_lamB.fa', 'w') as fh:
    for nm, sq in seqs[:400]: fh.write('>%s\n%s\n' % (nm, sq))
print('insertion-bearing reads:', len(seqs))
for k, v in sorted(junc.items(), key=lambda x: (x[0][1], x[0][2])):
    print('%-4s %-5s pos=%d  len~%d  reads=%d' % (k[0], k[1], k[2], k[3], v))
