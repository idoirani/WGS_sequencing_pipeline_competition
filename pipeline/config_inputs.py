"""Where THIS experiment's raw inputs live.

Separate from config.py on purpose: config.py describes the biological system
(reference, coordinates, thresholds) and rarely changes; this file describes one
experiment's files and changes every run.  00_sample_sheet.py is the only
consumer.
"""
EXP = ('/Users/ido/Dropbox/Bio/Phage/experiments/competition/'
       'CI857 deltalysis vs lysis competant 10212026')

# Lab plate record: one block per plate, rows = day, columns = condition.
CONDITIONS_XLSX = EXP + '/conditions_record.xlsx'

# Archive of basecalled libraries, named <run>_<n>_v<plate><row><col>.fastq.gz.
# Streamed member-by-member; never unpacked.
READ_ZIP = EXP + '/Sequencing_Days__2_3_4_rep_1_2.zip'

# Libraries sequenced in an earlier run that are not in the archive above.
# (plate, well_row, glob)  -- here: day 8 of plate 1, barcoded H1..H11.
EXTRA_RUNS = [('v1', 'H', EXP + '/H_sequncing_rep1/*_H*.fastq.gz')]
EXTRA_COL_RE = r'_H(\d+)\.fastq\.gz$'

# A library smaller than this is treated as a failed run (use=0).
MIN_FASTQ_BYTES = 5_000_000

# Which marker cassettes each CONDITION CLASS is built from.  This is a lab
# fact, written out so 17_concordance.py tests the data against the protocol
# rather than against itself.  First matching pattern wins (case-insensitive).
#   Ia  = aph(3')-Ia  kanR, lysis-competent   (TG5HVC_1)
#   IIa = aph(3')-IIa kanR, lysis-deficient   (DD7VVC_1)
#   cat = catA1       cmR,  lysis-competent   (DD7VVC_2)
CONDITION_MARKERS = [
    (r'100.*cmr',        ('cat',)),
    (r'100.*lysis\s*-',  ('IIa',)),
    (r'100',             ('Ia',)),
    (r'\+\s*/\s*-',      ('cat', 'IIa')),
    (r'\+\s*/\s*\+',     ('Ia', 'cat')),
    (r'double',          ('Ia', 'cat')),
]
