"""Single edit point for the whole pipeline.

Two things define a run:

  1. THIS FILE  - reference, coordinates, thresholds.  These are properties of
     the strains and the construct; they change only if you change the system.
  2. sample_sheet.csv - which libraries exist, what condition each one is, and
     where its reads live.  This changes every run.

Every numbered script imports from here; none of them contains a path, a sample
name or a coordinate.  Coordinates are 1-based inclusive on REF_FASTA.
"""
import os, csv, glob

# ------------------------------------------------------------------ 1. paths
ROOT      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORK      = os.path.join(ROOT, 'work')          # scratch: BAMs, probe FASTAs
OUT       = os.path.join(ROOT, 'handoff')       # per-stage JSON results
BAM       = os.path.join(WORK, 'bam', '{sample}.bam')

# Reference chromosome the reads are mapped to, and the contig to work on.
# This is the ancestor of the lysis-competent kanR strain; every coordinate
# below is on this assembly.
ANCESTOR_DIR = '/Users/ido/Dropbox/Bio/Strains/Chromosomes'
REF_STRAIN   = 'TG5HVC_1_lam_C1857_bor__kanR'
REF_FASTA    = os.path.join(ANCESTOR_DIR, REF_STRAIN, REF_STRAIN + '_1_reference.fasta')
CONTIG       = 'contig_1'

# The reads are NOT mapped to REF_FASTA. They are mapped to the four-variant
# reference built by 02_build_ref.py: REF_FASTA with the marker cassette masked
# out, plus one contig per marker allele. Anything that reads the BAMs -- most
# importantly the variant scan -- must use THIS file, or the contigs will not
# match the BAM header and bcftools will emit nothing at all.
MAP_FASTA    = os.path.join(WORK, 'ref_probe.fa')
# Four-variant locus reference (identical anchors, differing interiors) used by
# 04/05 for direct read assignment.  Built by 02b_build_locus_ref.py.
LOCUS_FASTA  = os.path.join(WORK, 'locus_ref.fa')

# The four parental assemblies.  Order matters only for reporting.
PARENTS = ['TG5HVC_1_lam_C1857_bor__kanR',         # kanR aph(3')-Ia,  lysis-competent
           'DD7VVC_2_WT+CI857_cmR',                # cmR  catA1,       lysis-competent
           'DD7VVC_1_cmS_kanR_CI857_delta_lysis',  # kanR aph(3')-IIa, lysis-deficient
           'TG5HVC_2_lam_WT']                      # lambda WT reference


# --------------------------------------------------------- 2. sample sheet
# sample_sheet.csv columns (extra columns are carried through untouched):
#   sample          library name; BAMs and outputs are keyed on it
#   fastq           reads.  Either a path/glob, or  zip:<archive>!<member>
#                   to stream a member straight out of a zip without unpacking
#   temp_C          competition temperature
#   day             serial-passage day (used by the time-series stages)
#   rep             biological replicate / plate
#   condition       free-text condition label from the lab's own record
#   group           condition group key, shared by the same condition across days
#   control_allele  non-empty for pure single-strain controls: the allele that
#                   library must show at fraction 1.0.  Drives the specificity
#                   check in 04 and the confusion matrix in 05.
#   lysis_plus_only 1 where the starting population was entirely lysis-competent
#   use             0 to exclude a library (failed run) without deleting the row
SHEET = os.environ.get('SAMPLE_SHEET', os.path.join(ROOT, 'sample_sheet.csv'))


def _truthy(v):
    return str(v).strip().lower() in ('1', 'y', 'yes', 'true')


def _load_sheet(path):
    if not os.path.exists(path):
        return []
    rows = []
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh):
            s = (r.get('sample') or '').strip()
            if not s or s.startswith('#'):
                continue
            if 'use' in r and not _truthy(r.get('use') or '1'):
                continue
            rows.append(r)
    return rows


SHEET_ROWS = _load_sheet(SHEET)
if not SHEET_ROWS:
    raise SystemExit('no usable rows in %s -- build the sample sheet first' % SHEET)

SAMPLES         = [r['sample'] for r in SHEET_ROWS]
ROW             = {r['sample']: r for r in SHEET_ROWS}
TEMP            = {r['sample']: int(float(r['temp_C'])) for r in SHEET_ROWS if r.get('temp_C')}
DAY             = {r['sample']: int(float(r['day'])) for r in SHEET_ROWS if r.get('day')}
REP             = {r['sample']: r.get('rep', '') for r in SHEET_ROWS}
CONDITION       = {r['sample']: r.get('condition', '') for r in SHEET_ROWS}
GROUP           = {r['sample']: r.get('group', '') for r in SHEET_ROWS}
# {sample: allele} for pure single-strain controls only
PURE_CONTROLS   = {r['sample']: r['control_allele'].strip()
                   for r in SHEET_ROWS if (r.get('control_allele') or '').strip()}
LYSIS_PLUS_ONLY = [r['sample'] for r in SHEET_ROWS if _truthy(r.get('lysis_plus_only') or '')]
FASTQ           = {r['sample']: r.get('fastq', '') for r in SHEET_ROWS}

# Controls grouped the other way round: allele -> [samples].  05 builds its
# confusion matrix from this, so adding control libraries needs no code change.
CONTROLS_BY_ALLELE = {}
for _s, _a in PURE_CONTROLS.items():
    CONTROLS_BY_ALLELE.setdefault(_a, []).append(_s)


# ------------------------------------------------------------- 3. prophage
# Span of the resident lambda prophage on CONTIG.
PROPHAGE = (3077400, 3123000)

# Two intervals used for the copy-number ratio in 06.  'lysis' is the endolysin
# R + holin S region; 'denom' is prophage sequence beside the cassette that is
# single-copy in the genome.  Both were checked unique by blastn before use.
REGION = {'lysis': (3106950, 3107733),
          'denom': (3108410, 3113337)}

# Repeat-like intervals excluded from the copy-number background.  Cryptic
# prophage / paralogue families; leaving them in biases the background median.
CRYPTIC = [(2246000, 2250000), (3349000, 3355000), (3365000, 3367000)]

# --------------------------------------------------------- 4. marker cassette
# The selectable-marker cassette sits at one locus; each parent carries a
# different allele there.  04 assigns reads to these four variants.
# name -> (strain, contig, cds_start, cds_end, label)
ALLELES = {
 'allele_Ia':  (PARENTS[0], CONTIG,     3105534, 3106349, "aph(3')-Ia  kanR, lysis+"),
 'allele_IIa': (PARENTS[2], 'contig_4',  430465,  431259, "aph(3')-IIa kanR, lysis-"),
 'allele_cat': (PARENTS[1], 'contig_5',  432913,  433572, "catA1       cmR,  lysis+"),
}
# Short allele keys as used inside 04/05 ('Ia', 'IIa', 'cat', 'bor').
VARIANTS = ['Ia', 'IIa', 'cat', 'bor']
# Cassette locus masked out of the backbone so reads are assigned by their
# interior sequence, not by shared flanks.
CASSETTE_MASK   = (3105000, 3106900)
CASSETTE_WINDOW = (3103500, 3109239)   # fetch window for locus reads
# Anchors that bracket the cassette.  Both are 1501 bp and byte-identical in
# every parent, so a read is assigned by its INTERIOR only.
ANCHOR_L = (3103500, 3105000)
ANCHOR_R = (3107739, 3109239)
FLANK    = 1500                        # decoy flank length in 02_build_ref.py

# ------------------------------------------- 5. cassette-independent lysis
# The lysis deletion removes endolysin R + holin S.  05 scores lysis status from
# this junction alone, so it does not depend on which marker a cell carries.
LYSIS_CONTEXT = (3106239, 3107738)     # R + S, present only in lysis+
LYSIS_ANCHOR  = (3107739, 3109239)     # shared right flank
LYSIS_MINUS   = ('DD7VVC_1_cmS_kanR_CI857_delta_lysis', 'contig_4', 430208, 431707)

# ----------------------------------------------------------- 6. host genes
# Receptor and its activator.  Strand matters for codon numbering in 09.
# Strand is taken from the strain's own annotation (both are minus here) and
# drives codon numbering in config.codon().  Getting it wrong silently reports
# every codon number mirrored about the middle of the gene.
GENES = {'lamB': {'span': (4324753, 4326093), 'strand': '-', 'codons': 446},
         'malT': {'span': (327940, 330645),  'strand': '-', 'codons': 901}}

# Tail-fibre locus, where prophage recombination with the Rac cryptic prophage
# shows up (12).
TFA = (3119000, 3120500)

# ------------------------------------------------------------- 7. thresholds
MAPQ        = 20      # primary-alignment MAPQ floor, all read counting
MIN_FLANK   = 200     # bp a read must carry either side of a junction
MIN_INS     = 200     # bp: smallest insertion treated as structural
MIN_DEL     = 50      # bp: smallest deletion treated as a lesion
MIN_JX_READS= 2       # reads needed to call an insertion junction
MIN_DISC    = 400     # bp inside the discriminating interior for locus assignment
BOOT        = 2000    # bootstrap replicates
# A pure control must be pure.  Fraction of its reads allowed to land on a
# DIFFERENT CASSETTE (cross-assignment onto wild-type is handled by the matrix
# itself and is not counted here).  A control above this is excluded from the
# confusion matrix and reported -- using a contaminated control as a matrix row
# biases every unmixed fraction in the run.
CTRL_MAX_OFF_CASSETTE = 0.10
MPILEUP_D   = 150     # bcftools mpileup per-file depth cap (unbounded stalls)
MIN_QUAL    = 30      # variant QUAL floor
MIN_DP      = 100     # joint-depth floor
THREADS     = int(os.environ.get('THREADS', '12'))

# Insertion-sequence families, called from the median inserted length.
# (name, expected bp, tolerance) -- confirmed by blastn of extracted inserts
# back to the strain's own chromosome.
IS_FAMILIES = [('IS1', 768, 40), ('IS5', 1195, 60), ('IS2', 1331, 70)]

os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.join(WORK, 'bam'), exist_ok=True)


def bam(sample):
    """Path to this sample's alignments.

    03_map.sh keeps the whole-genome BAM only when KEEP_BAM=1 (the variant scan
    needs it); otherwise it leaves the compact region BAM, which carries every
    interval stages 04-12 look at.  Stages do not care which one exists.
    """
    full = BAM.format(sample=sample)
    if os.path.exists(full):
        return full
    region = os.path.join(WORK, 'region', sample + '.bam')
    if os.path.exists(region):
        return region
    return full


def fastq_cmd(sample):
    """Shell command that writes this sample's reads to stdout.

    Supports a plain path, a glob, or 'zip:<archive>!<member>' so a 47 GB
    archive can be streamed member-by-member without ever being unpacked.
    """
    src = FASTQ[sample]
    if src.startswith('zip:'):
        arch, member = src[4:].split('!', 1)
        return 'unzip -p %s %s' % (_q(arch), _q(member))
    hits = sorted(glob.glob(src))
    if not hits:
        raise SystemExit('no fastq for %s at %s' % (sample, src))
    return 'cat ' + ' '.join(_q(h) for h in hits)


def _q(s):
    return "'" + s.replace("'", "'\\''") + "'"


def is_family(length):
    for name, centre, tol in IS_FAMILIES:
        if abs(length - centre) <= tol:
            return name
    return 'other'


def codon(gene, pos):
    """Chromosomal position -> codon number, honouring strand."""
    g = GENES[gene]
    lo, hi = g['span']
    return ((hi - pos) // 3 + 1) if g['strand'] == '-' else ((pos - lo) // 3 + 1)
