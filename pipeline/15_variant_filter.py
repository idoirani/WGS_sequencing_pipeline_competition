"""Turn the raw joint call set into the short list of real evolved variants.

Three filters, in order:
  1. parental    - sites where a non-reference ancestor already differs (stage 13)
  2. differential- a site must be >DIFF_HI in one sample and <DIFF_LO in another;
                   a variant shared by every sample is a reference/mapping
                   artefact, not something that evolved in this experiment
  3. repeat      - the +/-450 bp context must be single-copy in the reference.
                   LOAD-BEARING: the lambdoid cryptic prophages (tfa/tfaD,
                   packaging genes) form a 3-4 copy family at ~97% identity.
                   Reads from a diverged copy collapse onto the reference copy
                   and produce a dense block of 0.2-0.5 "variants" in exactly
                   those samples whose chromosome lacks the identical copy.
                   Without this filter that block dominates the result.
"""
import os, sys, json, subprocess, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from Bio import SeqIO

DIFF_HI, DIFF_LO, MIN_DP = 0.15, 0.03, 20
REP_PID, REP_LEN, FLANK = 85.0, 300, 450

vcf = os.path.join(C.WORK, 'par', 'evolved.vcf.gz')
q = subprocess.run(['bcftools', 'query', '-f', '%POS\t%REF\t%ALT[\t%AD]\n', vcf],
                   capture_output=True, text=True).stdout.strip().split('\n')
genome = {r.id: str(r.seq) for r in SeqIO.parse(C.MAP_FASTA, 'fasta')}[C.CONTIG]

sites = []
for line in q:
    f = line.split('\t')
    fr = []
    for ad in f[3:]:
        p = [int(x) if x != '.' else 0 for x in ad.split(',')]
        d = sum(p)
        fr.append(p[1] / d if d >= MIN_DP and len(p) > 1 else None)
    ok = [x for x in fr if x is not None]
    if ok and max(ok) > DIFF_HI and min(ok) < DIFF_LO:
        sites.append({'pos': int(f[0]), 'ref': f[1], 'alt': f[2], 'frac': fr})

# repeat test: one blast of all contexts at once
with tempfile.NamedTemporaryFile('w', suffix='.fa', delete=False) as fh:
    for s in sites:
        a, b = max(0, s['pos'] - FLANK - 1), s['pos'] + FLANK
        fh.write('>%d\n%s\n' % (s['pos'], genome[a:b]))
    probes = fh.name
hits = subprocess.run(['blastn', '-query', probes, '-subject', C.MAP_FASTA,
                       '-outfmt', '6 qseqid length pident'],
                      capture_output=True, text=True).stdout
ncopy = {}
for h in hits.strip().split('\n'):
    if not h:
        continue
    qid, ln, pid = h.split('\t')
    if int(ln) >= REP_LEN and float(pid) >= REP_PID:
        ncopy[int(qid)] = ncopy.get(int(qid), 0) + 1
for s in sites:
    s['n_copies'] = ncopy.get(s['pos'], 1)
    s['single_copy'] = s['n_copies'] <= 1

keep = [s for s in sites if s['single_copy']]
json.dump({'samples': C.SAMPLES, 'differential': sites, 'passing': keep},
          open(os.path.join(C.OUT, 'variants.json'), 'w'), indent=1)
print('differential %d | single-copy %d | repeat-masked %d'
      % (len(sites), len(keep), len(sites) - len(keep)))
for s in keep:
    print('%9d %-8s>%-8s  n=%d  %s' % (
        s['pos'], s['ref'][:8], s['alt'][:8], s['n_copies'],
        ' '.join('%4s' % ('%.2f' % x if x is not None else '  . ') for x in s['frac'])))

# ---------------------------------------------------------------------------
# classify what survives.  Two nuisance classes remain after repeat masking:
#   homopolymer - indel inside a >=5 bp mononucleotide run.  ONT's dominant
#                 error mode; appears at low fraction in many samples at once,
#                 including samples that cannot carry the mutation.
#   ubiquitous  - non-zero in most samples including the pure single-strain
#                 controls, so it did not arise during the experiment
#                 (reference base error, or a paralog below the repeat cutoff).
# What is left is sample-specific: high in one or two samples, absent in the
# rest.  Those are the real evolved variants.
# ---------------------------------------------------------------------------
import re, csv

def homopolymer(seq, pos, k=5):
    w = seq[pos - 1 - 12: pos - 1 + 12]
    return bool(re.search(r'(A{%d,}|C{%d,}|G{%d,}|T{%d,})' % ((k,) * 4), w))

rows = []
for s in keep:
    fr = [x for x in s['frac'] if x is not None]
    nz = sum(1 for x in fr if x > 0.05)
    hp = homopolymer(genome, s['pos']) and len(s['ref']) != len(s['alt'])
    hi = max(fr)
    rest = sorted(fr)[-2] if len(fr) > 1 else 0.0     # 2nd highest sample
    # Sample-specificity is the discriminator and it OUTRANKS homopolymer
    # context: the ONT indel floor is low fraction spread over many samples
    # (measured here at 4.9%, 19/392 reads in the 30 C controls), whereas a
    # real indel is high in one sample and absent from the rest.  A genuine
    # frameshift inside a T-run looks nothing like the floor it sits on.
    specific = hi > DIFF_HI and rest < 0.05
    if specific and (hi >= 0.30 or not hp):
        cls = 'EVOLVED'
    elif hp:
        cls = 'homopolymer (ONT error floor)'
    elif nz > len(fr) / 2:
        cls = 'ubiquitous (not evolved)'
    else:
        cls = 'low-confidence'
    s['klass'] = cls
    rows.append(s)

hdr = ['pos', 'ref', 'alt', 'n_copies', 'class'] + C.SAMPLES
with open(os.path.join(C.WORK, 'variants_genomewide.csv'), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(hdr)
    for s in rows:
        w.writerow([s['pos'], s['ref'], s['alt'], s['n_copies'], s['klass']] +
                   ['' if x is None else round(x, 4) for x in s['frac']])
json.dump({'samples': C.SAMPLES, 'differential': sites, 'passing': rows},
          open(os.path.join(C.OUT, 'variants.json'), 'w'), indent=1)

for cls in ('EVOLVED', 'homopolymer (ONT error floor)', 'ubiquitous (not evolved)'):
    sel = [s for s in rows if s['klass'] == cls]
    print('\n%-30s n=%d' % (cls, len(sel)))
    for s in sel:
        print('  %9d %-9s>%-9s %s' % (s['pos'], s['ref'][:9], s['alt'][:9],
              ' '.join('%4s' % ('%.2f' % x if x is not None else '  . ') for x in s['frac'])))
