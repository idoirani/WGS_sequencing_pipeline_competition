#!/usr/bin/env python
"""Gene-wide loss-of-function calling for the receptor and its activator.

Stages 09/11 answered a day-8 question and only knew about the lesions that had
already been found there, at hard-coded positions.  This stage scans the WHOLE
coding sequence of each gene in config.GENES and classifies every spanning read
into exactly one class, so a new experiment finds new lesions without edits:

  IS     >= MIN_INS bp insertion inside the CDS      (IS family by length)
  DEL    >= MIN_DEL bp deletion inside the CDS
  BREAK  >= MIN_INS bp clip starting inside the CDS  (read diverges there)
  FS     short indel at a site that passes the error-floor test -> frameshift
  NS     substitution at such a site that creates a stop codon
  WT     spans the whole CDS with none of the above

One read gets one class, and every class shares one denominator, so the
fractions are comparable and sum to 1.

Point lesions are NOT trusted from a single read.  A candidate site is kept
only if some library carries it at >= POINT_MIN_F while the pure 30 C control
libraries -- which cannot have evolved it -- stay below POINT_FLOOR_MULT times
that library's level.  This is what separates a real frameshift from the ONT
homopolymer indel error class.  Crucially the sample-specificity test runs
BEFORE any homopolymer test: a genuine frameshift often sits in a mononucleotide
run, and filtering on context first discards it.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
import subprocess, re, json, math, csv, collections
from Bio import SeqIO
from Bio.Seq import Seq

CIG = re.compile(r'(\d+)([MIDNSHP=X])')
S = C.SAMPLES
PAD = 2000            # bp either side of the CDS to fetch reads from
FLANK_REQ = 150       # bp a read must carry outside the CDS to count as spanning
POINT_MIN_F = 0.10    # a candidate point site must reach this in some library
# ...and the libraries that CANNOT have evolved it (pure, 30 C) must be near
# silent there.  A purely relative test is not enough: a 1 bp deletion sitting
# at 5.8% in every control still cleared "floor < half the peak" and was then
# counted as a frameshift in all 65 libraries, putting 9-30% false disruption
# into pure controls.  The cap is per event kind because the ONT error classes
# differ by an order of magnitude -- homopolymer indels dominate, substitutions
# are rarer -- and the peak must also stand well clear of whatever floor there is.
POINT_MAX_CONTROL = {'sub': 0.02, 'ins': 0.03, 'del': 0.03}
POINT_PEAK_OVER_FLOOR = 4.0
# A fraction alone is not evidence.  In a 4-read library one sequencing error
# is f = 0.25 and clears any fraction threshold, which fills the catalogue with
# edge artefacts at the ends of thinly covered genes.  The peak library must
# also carry real read support.
POINT_MIN_READS = 5   # supporting reads in the library where the site peaks
POINT_MIN_DEPTH = 20  # ...out of at least this many reads at the gene
# Only loss-of-function is reported: f_disrupted counts nonsense and frameshift.
# Missense calls are not actionable here and would swamp the catalogue.
KEEP_EFFECTS = ('nonsense', 'frameshift')


def wilson(x, n, z=1.96):
    if not n:
        return (None, None)
    d = 1 + z * z / n
    c = (x / n + z * z / (2 * n)) / d
    h = z * math.sqrt((x / n) * (1 - x / n) / n + z * z / (4 * n * n)) / d
    return max(0., c - h), min(1., c + h)


genome = {r.id: str(r.seq) for r in SeqIO.parse(C.REF_FASTA, 'fasta')}
CHROM = genome[C.CONTIG]

# control libraries that cannot have evolved a 37 C lesion: pure, 30 C
FLOOR_SAMPLES = [s for s in C.PURE_CONTROLS if C.TEMP.get(s) == 30] or list(C.PURE_CONTROLS)


def read_events(sample, gene):
    """-> (list of per-read records, n_fetched).  No classification yet."""
    lo, hi = C.GENES[gene]['span']
    p = subprocess.run(['samtools', 'view', '-F', '0x904', '-q', str(C.MAPQ),
                        C.bam(sample), '%s:%d-%d' % (C.CONTIG, lo - PAD, hi + PAD)],
                       capture_output=True, text=True).stdout
    out = []
    for line in p.split('\n'):
        if not line:
            continue
        f = line.split('\t', 10)
        pos, cig, seq = int(f[3]), f[5], f[9]
        if cig == '*':
            continue
        r, q = pos, 0
        big_ins, big_del, clips, pts = [], [], [], []
        for i, (n, o) in enumerate(CIG.findall(cig)):
            n = int(n)
            if o == 'I':
                if n >= C.MIN_INS and lo <= r <= hi:
                    big_ins.append((r, n))
                elif 1 <= n < C.MIN_INS and lo <= r <= hi:
                    pts.append((r, 'ins', n))
                q += n
            elif o in 'M=X':
                # substitutions against the reference
                for k in range(n):
                    rp = r + k
                    if lo <= rp <= hi:
                        b = seq[q + k:q + k + 1].upper()
                        if b and b != CHROM[rp - 1].upper():
                            pts.append((rp, 'sub', b))
                r += n; q += n
            elif o in 'DN':
                if n >= C.MIN_DEL and lo <= r <= hi:
                    big_del.append((r, n))
                elif 1 <= n < C.MIN_DEL and lo <= r <= hi:
                    pts.append((r, 'del', n))
                r += n
            elif o in 'SH':
                if n >= C.MIN_INS:
                    clips.append(pos if i == 0 else r)
                if o == 'S':
                    q += n
        end = r
        out.append({'start': pos, 'end': end, 'ins': big_ins, 'del': big_del,
                    'clip': [c for c in clips if lo <= c <= hi],
                    'pts': pts,
                    'spans': pos <= lo - FLANK_REQ and end >= hi + FLANK_REQ,
                    'enters': pos <= lo - FLANK_REQ or end >= hi + FLANK_REQ})
    return out


def annotate_point(gene, pos, kind, val):
    """-> ('nonsense'|'frameshift'|'missense'|'inframe'|'silent', codon_no)."""
    g = C.GENES[gene]
    cn = C.codon(gene, pos)
    if kind in ('ins', 'del'):
        return ('frameshift' if val % 3 else 'inframe'), cn
    lo, hi = g['span']
    cds = CHROM[lo - 1:hi]
    if g['strand'] == '-':
        off = hi - pos
        alt = str(Seq(val).complement())
    else:
        off = pos - lo
        alt = val
    seqc = cds if g['strand'] == '+' else str(Seq(cds).reverse_complement())
    i = off
    if i >= len(seqc):
        return 'missense', cn
    mut = seqc[:i] + alt + seqc[i + 1:]
    ci = (i // 3) * 3
    old_aa = str(Seq(seqc[ci:ci + 3]).translate())
    new_aa = str(Seq(mut[ci:ci + 3]).translate())
    if new_aa == '*' and old_aa != '*':
        return 'nonsense', cn
    return ('silent' if new_aa == old_aa else 'missense'), cn


def main():
    # ---- pass 1: collect every read's raw events ---------------------------
    data = {}
    for gene in C.GENES:
        for s in S:
            data[(s, gene)] = read_events(s, gene)

    # ---- decide which point sites are real --------------------------------
    # A read contributes point candidates only if it is also in the denominator.
    # Counting events from partially-overlapping reads against a denominator of
    # spanning reads produced fractions above 1.
    def counted(r):
        return r['spans'] or r['ins'] or r['del'] or r['clip']

    cand = collections.defaultdict(dict)   # (gene,pos,kind,val) -> {sample: (x, n)}
    LEN = collections.defaultdict(collections.Counter)   # modal indel length
    denom = {}
    for (s, gene), reads in data.items():
        use = [r for r in reads if counted(r)]
        n = len(use)
        denom[(s, gene)] = n
        if not n:
            continue
        cnt = collections.Counter()
        for r in use:
            # An indel is keyed by position and kind only, NOT by its exact
            # length: ONT reports one true 25 bp insertion as 24-27 bp, and
            # keying on length splits a single real lesion into four sites,
            # each below threshold.  Substitutions keep their base.
            ev = {(pos, kind, (val if kind == 'sub' else '*')) for pos, kind, val in r['pts']}
            for pos, kind, val in ev:
                cnt[(gene, pos, kind, val)] += 1
            for pos, kind, val in r['pts']:
                if kind != 'sub':
                    LEN[(gene, pos, kind)][val] += 1
        for k, c in cnt.items():
            cand[k][s] = (c, n)

    real, borderline = {}, []
    for k, per in cand.items():
        peak_s = max(per, key=lambda s: per[s][0] / per[s][1])
        x, n = per[peak_s]
        peak = x / n
        if peak < POINT_MIN_F or x < POINT_MIN_READS or n < POINT_MIN_DEPTH:
            continue
        gene, pos, kind, val = k
        size = (LEN[(gene, pos, kind)].most_common(1)[0][0]
                if kind != 'sub' and LEN[(gene, pos, kind)] else val)
        eff, cn = annotate_point(gene, pos, kind, size)
        if eff not in KEEP_EFFECTS:
            continue
        floor = max([per[c][0] / per[c][1] for c in FLOOR_SAMPLES if c in per] or [0.0])
        # Sample-specificity FIRST: a real lesion is high in one library and
        # near-absent in libraries that cannot carry it.  Running a homopolymer
        # test before this one discards genuine frameshifts in mononucleotide
        # runs -- which is how the day-8 receptor frameshift was nearly lost.
        rec = {'effect': eff, 'codon': cn, 'peak_f': peak, 'floor_f': floor,
               'peak_sample': peak_s, 'peak_x': x, 'peak_n': n}
        if floor > POINT_MAX_CONTROL.get(kind, 0.02) or peak < POINT_PEAK_OVER_FLOOR * floor:
            borderline.append((k, rec))
            continue
        rec['size'] = size
        real[k] = rec
    borderline.sort(key=lambda kv: -kv[1]['peak_f'])

    # ---- pass 2: one class per read ---------------------------------------
    rows, cat = [], []
    for k, v in sorted(real.items()):
        gene, pos, kind, val = k
        cat.append({'gene': gene, 'position': pos, 'kind': kind, 'change': (val if kind == 'sub' else '%d bp' % v['size']),
                    'effect': v['effect'], 'codon': v['codon'],
                    'peak_sample': v['peak_sample'],
                    'peak_reads': '%d/%d' % (v['peak_x'], v['peak_n']),
                    'peak_fraction': round(v['peak_f'], 4),
                    'control_fraction': round(v['floor_f'], 4)})

    for gene in C.GENES:
        for s in S:
            reads = data[(s, gene)]
            cnt = collections.Counter()
            junc = collections.Counter()
            for r in reads:
                if r['ins']:
                    cnt['IS'] += 1
                    p, L = min(r['ins'])
                    junc[(p // 60 * 60, C.is_family(L))] += 1
                elif r['del']:
                    cnt['DEL'] += 1
                    junc[(min(r['del'])[0] // 60 * 60, 'del')] += 1
                elif r['clip'] and r['enters']:
                    cnt['BREAK'] += 1
                    junc[(min(r['clip']) // 60 * 60, 'break')] += 1
                else:
                    ev = {(p, kd, (vl if kd == 'sub' else '*')) for p, kd, vl in r['pts']}
                    hit = [real[(gene, p, kd, vl)]['effect']
                           for p, kd, vl in ev if (gene, p, kd, vl) in real]
                    if 'nonsense' in hit:
                        cnt['NS'] += 1
                    elif 'frameshift' in hit:
                        cnt['FS'] += 1
                    elif r['spans']:
                        cnt['WT'] += 1
            n = sum(cnt.values())
            dis = n - cnt['WT']
            lo, hi = wilson(dis, n)
            rows.append({'sample': s, 'gene': gene, 'day': C.DAY.get(s, ''),
                         'rep': C.REP.get(s, ''), 'temp_C': C.TEMP.get(s, ''),
                         'group': C.GROUP.get(s, ''), 'condition': C.CONDITION.get(s, ''),
                         'n_reads': n, 'n_IS': cnt['IS'], 'n_DEL': cnt['DEL'],
                         'n_BREAK': cnt['BREAK'], 'n_FS': cnt['FS'], 'n_NS': cnt['NS'],
                         'n_WT': cnt['WT'],
                         'f_disrupted': round(dis / n, 4) if n else '',
                         'ci_lo': round(lo, 4) if n else '', 'ci_hi': round(hi, 4) if n else '',
                         'f_IS': round(cnt['IS'] / n, 4) if n else '',
                         'f_DEL': round(cnt['DEL'] / n, 4) if n else '',
                         'f_BREAK': round(cnt['BREAK'] / n, 4) if n else '',
                         'f_point': round((cnt['FS'] + cnt['NS']) / n, 4) if n else '',
                         'n_junctions': len(junc)})

    out1 = os.path.join(C.ROOT, 'lesions_per_sample.csv')
    with open(out1, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    out2 = os.path.join(C.ROOT, 'lesion_catalogue.csv')
    with open(out2, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(cat[0]) if cat else
                           ['gene', 'position', 'kind', 'change', 'effect', 'codon',
                            'peak_sample', 'peak_reads', 'peak_fraction',
                            'control_fraction'])
        w.writeheader(); w.writerows(cat)
    json.dump({'rows': rows, 'catalogue': cat}, open(os.path.join(C.OUT, 'lesions.json'), 'w'))

    print('point-lesion catalogue (%d sites kept of %d candidates):' % (len(cat), len(cand)))
    for c in cat:
        print('  %-5s %-9d %-4s %-6s %-11s codon %-5d  %-7s %-8s peak %.3f control %.3f'
              % (c['gene'], c['position'], c['kind'], str(c['change'])[:6],
                 c['effect'], c['codon'], c['peak_sample'], c['peak_reads'],
                 c['peak_fraction'], c['control_fraction']))
    if borderline:
        print('\nrejected as error floor (reported, not counted as disruption):')
        for (gene, pos, kind, val), v in borderline[:10]:
            print('  %-5s %-9d %-4s %-6s %-11s codon %-5d  %-7s peak %.3f control %.3f'
                  % (gene, pos, kind, str(val)[:6], v['effect'], v['codon'],
                     v['peak_sample'], v['peak_f'], v['floor_f']))
    print('\nwrote %s, %s' % (out1, out2))


if __name__ == '__main__':
    main()
