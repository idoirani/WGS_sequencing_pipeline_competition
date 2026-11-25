#!/usr/bin/env python
"""Build the mapping reference.

Backbone : REF_STRAIN closed assembly, cassette locus hard-masked to N.
Decoys   : the three cassette alleles, each with FLANK bp of its own native
           flank, so every allele is measured in an equivalent sequence context.

Writes work/ref_probe.fa (= config.MAP_FASTA) and work/ref_probe.json.
Everything that later reads a BAM must use THIS fasta: the BAM header carries
the decoy contigs, and pointing a variant caller at the plain ancestor instead
produces no output at all because the contig names do not match.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from Bio import SeqIO
from Bio.Seq import Seq


def load(strain):
    p = os.path.join(C.ANCESTOR_DIR, strain, strain + '_1_reference.fasta')
    return {r.id: r for r in SeqIO.parse(p, 'fasta')}


def main():
    recs, meta = [], {}

    bb = load(C.REF_STRAIN)
    s, e = C.CASSETTE_MASK
    seq = str(bb[C.CONTIG].seq)
    bb[C.CONTIG].seq = Seq(seq[:s] + 'N' * (e - s) + seq[e:])
    for r in bb.values():
        r.description = ''
        recs.append(r)
    meta['backbone'] = {'strain': C.REF_STRAIN,
                        'masked': [C.CONTIG, s, e],
                        'contigs': {k: len(v.seq) for k, v in bb.items()}}

    for name, (strain, contig, cs, ce, gene) in C.ALLELES.items():
        src = load(strain)[contig]
        sub = src.seq[cs - 1 - C.FLANK: ce + C.FLANK]
        recs.append(SeqIO.SeqRecord(sub, id=name, description=''))
        meta[name] = {'gene': gene, 'strain': strain, 'source_contig': contig,
                      'source_span': [cs, ce], 'length': len(sub),
                      # CDS interior in decoy coordinates, 1-based inclusive
                      'cds': [C.FLANK + 1, C.FLANK + (ce - cs + 1)]}

    os.makedirs(C.WORK, exist_ok=True)
    SeqIO.write(recs, C.MAP_FASTA, 'fasta')
    json.dump(meta, open(os.path.join(C.WORK, 'ref_probe.json'), 'w'), indent=1)
    print('contigs:', len(recs), 'total bp:', sum(len(r.seq) for r in recs))
    for k in C.ALLELES:
        print(' ', k, meta[k]['gene'], 'len', meta[k]['length'], 'cds', meta[k]['cds'])


if __name__ == '__main__':
    main()
