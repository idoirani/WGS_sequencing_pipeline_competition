#!/usr/bin/env python
"""Build the four-variant LOCUS reference used by 04 and 05.

Each variant is  ANCHOR_L + <that parent's interior> + ANCHOR_R, where both
anchors are taken from the reference strain and are therefore byte-identical
across all four entries.  A read can only be assigned by its INTERIOR, which is
the whole point: with differing flanks, assignment would partly reflect flank
similarity rather than which cassette the cell actually carries.

The interior of each parent is located by blasting the two anchors against that
parent's own assembly, so a parent whose locus sits on a different contig, or in
the opposite orientation, is handled without hand-entered coordinates.

Writes work/locus_ref.fa (= config.LOCUS_FASTA) and work/locus_ref.json.
"""
import os, sys, json, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from Bio import SeqIO
from Bio.Seq import Seq

P = {'Ia': C.PARENTS[0], 'cat': C.PARENTS[1], 'IIa': C.PARENTS[2], 'bor': C.PARENTS[3]}
MIN_ANCHOR_LEN = 1490
MIN_ANCHOR_PID = 99.5


def genome(strain):
    p = os.path.join(C.ANCESTOR_DIR, strain, strain + '_1_reference.fasta')
    return {r.id: str(r.seq) for r in SeqIO.parse(p, 'fasta')}


def main():
    os.makedirs(C.WORK, exist_ok=True)
    G = {k: genome(v) for k, v in P.items()}

    al0, al1 = C.ANCHOR_L
    ar0, ar1 = C.ANCHOR_R
    AL = G['Ia'][C.CONTIG][al0 - 1:al1]
    AR = G['Ia'][C.CONTIG][ar0 - 1:ar1]
    anchors = os.path.join(C.WORK, 'anchors.fa')
    open(anchors, 'w').write('>AL\n%s\n>AR\n%s\n' % (AL, AR))

    rec = {}
    for k, strain in P.items():
        tmp = os.path.join(C.WORK, '_g_%s.fa' % k)
        open(tmp, 'w').write(''.join('>%s\n%s\n' % (c, s) for c, s in G[k].items()))
        hit = {}
        out = subprocess.run(
            ['blastn', '-query', anchors, '-subject', tmp, '-outfmt',
             '6 qseqid sseqid sstart send length pident', '-evalue', '1e-30'],
            capture_output=True, text=True, check=True).stdout
        for ln in out.splitlines():
            q, sc, ss, se, L, pid = ln.split('\t')
            if int(L) >= MIN_ANCHOR_LEN and float(pid) > MIN_ANCHOR_PID and q not in hit:
                hit[q] = (sc, int(ss), int(se))
        os.remove(tmp)
        if 'AL' not in hit or 'AR' not in hit:
            raise SystemExit('%s: anchors not found (AL=%s AR=%s)'
                             % (k, 'AL' in hit, 'AR' in hit))
        (cl, l1, l2), (cr, r1, r2) = hit['AL'], hit['AR']
        if cl != cr:
            raise SystemExit('%s: anchors on different contigs (%s, %s)' % (k, cl, cr))
        fwd = l1 < l2
        if fwd:
            interior = G[k][cl][l2:r1 - 1]
        else:
            interior = str(Seq(G[k][cl][r2:l1 - 1]).reverse_complement())
        rec[k] = {'strain': strain, 'contig': cl, 'anchorL': [l1, l2],
                  'anchorR': [r1, r2], 'fwd': fwd, 'interior_len': len(interior),
                  'locus_len': len(AL) + len(interior) + len(AR), 'seq': AL + interior + AR}
        print('%-4s %-9s anchorL %d-%d anchorR %d-%d fwd=%s interior=%d bp'
              % (k, cl, l1, l2, r1, r2, fwd, len(interior)))

    with open(C.LOCUS_FASTA, 'w') as fh:
        for k in C.VARIANTS:
            fh.write('>var_%s\n%s\n' % (k, rec[k]['seq']))
    json.dump({k: {x: y for x, y in v.items() if x != 'seq'} for k, v in rec.items()},
              open(os.path.join(C.WORK, 'locus_ref.json'), 'w'), indent=1)
    print('\nlocus_ref.fa:', {k: rec[k]['locus_len'] for k in rec},
          ' anchors', len(AL), '/', len(AR))


if __name__ == '__main__':
    main()
