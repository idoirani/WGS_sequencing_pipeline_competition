# Serial-passage time course of a λ *cI857* lysogen competition

**65 Oxford Nanopore whole-genome libraries, 11 condition groups, passage days 1–8,
two replicate plates.**

---

## 1. Executive summary

1. **Induction drives a 10–50-fold burst of phage DNA that then collapses.** At
   37 °C, phage genome equivalents per chromosome peak at day 1–2 (group means
   4.1–26.3; the single highest library is 49.7) and fall to 0.8–2.9 by day 8.
   At 30 °C the same quantity never leaves 0.97–2.9 at any time point. This is
   the clearest and most reproducible signal in the dataset.

2. **Receptor loss is late, large, and strictly temperature-dependent.** The
   *lamB* disrupted fraction is at the control floor through day 2 (37 °C group
   means 0.000–0.036), rises to 0.34–0.46 at day 4, and reaches 0.76–0.87 at day
   8 — in **every** 37 °C group, including those with no lysis-deficient
   competitor at all. At 30 °C it never exceeds 0.044 at any time point.

3. **The maltose activator *malT* is hit too, and earlier, but transiently.**
   37 °C group means rise from 0.000–0.018 at day 1 to 0.30–0.57 at day 4, then
   fall to 0.07–0.31 at day 8 while *lamB* keeps climbing. A single 25 bp
   insertion at codon 341 accounts for most of it (§5.2).

4. **Neither competitor is eliminated.** In the ⁺:⁻ mixtures the lysis-deficient
   fraction stays within roughly ±0.2 of its inoculum at every time point in both
   temperatures, with no monotone trend. The cmR competitor does gain in the ⁺:⁺
   mixtures by day 8 (0.79 and 0.76 cmR), but those are single flasks.

5. **Seven libraries disagree with their condition label**, including two that
   look like a transposed pair and one pure control that is 18.5 % contaminated.
   All are reported, none relabelled (§6).

---

## 2. Design and data

| | |
|---|---|
| Libraries | 66 sequenced, 65 analysed (v1B1 excluded, 191 kB) |
| Plates | v1 = replicate 1, v2 = replicate 2 |
| Days | rep 1: 2, 4, 8 · rep 2: 1, 2, 4 |
| Conditions | 11 groups: 3 pure controls, ⁺:⁻ and ⁺:⁺ mixtures at 50/50 and 20/80, double lysogens, each at 30 °C and 37 °C |
| Depth | median chromosome depth 39–355×, 100 % of reads mapped in every library |

Replication is the binding limitation. Days 2 and 4 have both plates; day 1 has
only replicate 2 and day 8 only replicate 1. **No group-day cell has more than
two flasks**, so nothing here is tested for significance; every figure shows each
flask individually and every aggregated row carries `n_flasks` and the min/max
across flasks rather than a standard error.

## 3. Method

Reads were streamed directly from the 47 GB archive (never unpacked) and mapped
with `minimap2 -ax map-ont --secondary=no` onto a reference consisting of the
ancestral chromosome with the marker cassette masked, plus one decoy contig per
cassette allele. Each library was then reduced to a region BAM and a 1 kb depth
table; whole-genome BAMs were discarded.

**Nothing is called from an assembly.** A consensus assembly of a mixed
population reports the majority allele as if it were the only allele — which is
exactly the quantity being measured here.

Marker fractions come from direct read assignment at a four-variant locus whose
flanking anchors are byte-identical across all four parents, so a read is
assigned on cassette sequence alone. Lysis genotype is read independently at the
*S*–*R*–*Rz* junction. Copy number is the median prophage depth over the
genome-wide background, with a faux-prophage null for the error bar.

### 3.1 Cross-assignment is calibrated per plate row, not pooled

Near-identical alleles cross-assign, so pure single-strain controls are used as
the rows of a confusion matrix which is inverted analytically. The measured
IIa→wild-type leak is **not** a constant of the assay: it ranges 0.055–0.139
between plate rows. Pooling all controls would push that spread into every
mixture as a systematic error, so each library is unmixed with the controls from
its own replicate and day (the plate layout places the three pure controls in
columns 1–3 of every row). Where a row's control is missing or failed QC, the
pooled mean is substituted and the substitution is recorded.

Matching tightened the allele budget — the sum of the three unmixed fractions,
which should be 1.0 — from **0.944–1.091 to 0.963–1.077** across all 65
libraries, while leaving the day-8 values unchanged.

### 3.2 Validation against the day-8 analysis

The pipeline was re-run end to end on the 11 previously analysed day-8 libraries
before any new library was trusted. **All 121 checked metrics reproduce the saved
day-8 tables exactly** (`regression_day8.csv`): read counts, marker fractions,
allele budgets, lysis-junction fractions, copy number, and *lamB*/*malT*
disruption. The only deviations are in the fourth decimal and come from rounding
in the stored CSV.

---

## 4. Results

### 4.1 Prophage copy number (Figure 7)

| Group | 30 °C d1→d8 | | Group | 37 °C d1→d8 |
|---|---|---|---|---|
| 50/50 ⁺:⁻ | 1.61 → 1.33 | | 50/50 ⁺:⁻ | 12.93 → **26.30** (d2) → 9.00 → 2.87 |
| 20/80 ⁺:⁻ | 1.19 → 1.18 | | 20/80 ⁺:⁻ | 10.34 → 13.57 → 8.66 → 2.01 |
| double lysogen | 2.94 → 1.51 | | 50/50 ⁺:⁺ | 6.60 → 9.54 → 2.42 → 1.93 |
| | | | 20/80 ⁺:⁺ | 4.12 → **18.04** (d2) → 2.49 → 2.02 |
| | | | double lysogen | 7.59 → 12.02 → 7.56 → 0.76 |

At 30 °C the pure single lysogens sit at 0.97–1.09 — one integrated copy, as they
should — and the double lysogens at 1.51–2.94, consistent with two. At 37 °C
every group rises an order of magnitude and then falls back toward one by day 8.

The measured quantity is phage genome equivalents per chromosome, which after
induction counts replicating phage DNA as well as integrated prophage. The day
1–2 peak is therefore the induction burst itself, not a change in integrated
copy number, and the decline by day 8 is the population returning to a
lysogen-dominated state.

### 4.2 Marker and lysis composition (Figure 8)

In the ⁺:⁻ mixtures the lysis-deficient fraction oscillates around its inoculum
without resolving: 50/50 at 37 °C runs 0.46 → 0.31 → 0.64 → 0.57, and 20/80 runs
0.77 → 0.41 → 0.79 → 0.51. At 30 °C the same groups run 0.35 → 0.39 → 0.53 → 0.59
and 0.74 → 0.54 → 0.74 → 0.78. The day-to-day swings are as large as the
difference between the two inocula, which with one or two flasks per point means
**no directional conclusion is supportable** — including the day-8-only reading
reported previously.

The double lysogens hold lysis-deficient at 0.000–0.010 throughout, as they must,
and the pure controls return 0.000 and 1.000 exactly. These are the internal
checks that the junction assay is doing what it claims.

In the ⁺:⁺ mixtures — where both competitors are lysis-competent and differ only
in marker — the cmR strain reaches 0.79 (50/50) and 0.76 (20/80) by day 8 from
starting points near 0.5 and 0.37. Both are single flasks with no day-8 replicate.

### 4.3 Receptor and activator loss (Figure 9)

Mean disrupted fraction across all 37 °C libraries:

| Gene | day 1 (n=5) | day 2 (n=10) | day 4 (n=10) | day 8 (n=5) |
|---|---|---|---|---|
| *lamB* | 0.015 | 0.021 | 0.391 | **0.826** |
| *malT* | 0.004 | 0.146 | **0.452** | 0.268 |

*lamB* rises monotonically and is still rising at day 8. *malT* peaks at day 4
and then declines — though day 8 has half as many libraries as day 4, so the
decline rests on five flasks.

Receptor loss reaches 0.76–0.87 in **all five** day-8 37 °C libraries, including
the two with no lysis-deficient competitor and the double lysogen. Whatever is
driving receptor loss is a consequence of induction, not of competing against a
lysis-deficient partner.

---

## 5. What the lesions are

### 5.1 Mechanism is mostly structural, not point mutation

Across the libraries above 20 % disruption, *lamB* is broken predominantly by
insertion-sequence insertions and their junction breakpoints. Point lesions
matter in only two libraries: the Trp99→stop in v1H7 and the codon-115
frameshift in v1H11, both of which were already known from the day-8 analysis and
are recovered here independently.

One new *lamB* nonsense mutation appears in replicate 2: **codon 266 in v2D8**,
at 0.238 and absent from every control.

### 5.2 A convergent 25 bp insertion in *malT*

A 25 bp insertion at *malT* codon 341 appears in **12 of 65 libraries at >5 %**,
reaching 0.284. It is present **only** in 37 °C libraries, in both replicates, in
multiple condition groups, and is at exactly **0.000 in all 17 pure 30 °C control
libraries**. The observed insert length varies 24–27 bp, which is ONT length
jitter around one true event — the caller originally split it into four
sub-threshold sites and was corrected to key indels by position rather than
length.

Independent recurrence in both replicates and across condition groups makes this
convergent evolution at the maltose regulon rather than a single lineage
spreading. *malT* activates *lamB* transcription, so losing *malT* removes the
receptor without touching the receptor gene — a second, cheaper route to the same
phenotype, and one that appears earlier than *lamB* disruption itself.

### 5.3 Two candidate lesions were rejected

A *malT* nonsense at codon 524 (peak 0.152, control 0.030) and a *malT* 1 bp
frameshift at codon 619 (peak 0.135, control 0.058) both fail the control-floor
test and are reported as unresolved rather than counted. Before that filter was
tightened, the 1 bp deletion alone put 9–43 % false disruption into pure 30 °C
control libraries. The residual false-disruption floor is now 0.000–0.026 in
every control except one 18-read library at 0.111.

---

## 6. Libraries that disagree with their label

58 of 65 libraries match the condition they are labelled with. The seven that do
not are listed in `genotype_concordance.csv` and reported here **without
relabelling** — a disagreement can be a plating error, a pipetting error, or a
transposed barcode, and only the lab can say which.

| Library | Labelled | Observed | Reading |
|---|---|---|---|
| **v2D1** | 100 % kanR lysis⁺ control | 80.1 % Ia, **16.7 % cat** | Contaminated control. Excluded from every confusion matrix. |
| **v2B8** | 20/80 ⁺:⁺ 37 °C | IIa+cat, lysis⁺ 0.735 | Behaves as a ⁺:⁻ flask |
| **v2B9** | 50/50 ⁺:⁻ 37 °C | Ia+cat, lysis⁺ 0.9997 | Behaves as a ⁺:⁺ flask — v2B8/v2B9 look like a transposed pair |
| **v2D10** | 20/80 ⁺:⁻ 37 °C | IIa 1.02, cat 0.031, lysis⁺ 0.049 | Probably biological: the cmR competitor was lost |
| v1B5, v1B9, v1H5 | ⁺:⁻ mixtures | trace Ia at 0.09–0.15 | Marginal; all are low-coverage (43–172 locus reads) |

**The day-8 transposition is independently confirmed.** The lab's own plate
record already lists columns 8 and 10 swapped in the day-8 row relative to every
other row, and the read data agree: v1H8 behaves as a ⁺:⁻ flask and v1H10 as ⁺:⁺,
exactly as the record states. The record was used as written.

---

## 7. Limitations

- **n ≤ 2 everywhere.** No statistical test is reported and none would be
  meaningful. Day 1 is replicate 2 only; day 8 is replicate 1 only, so the
  day-1-to-day-8 comparison is also a between-replicate comparison.
- Copy number after induction mixes integrated prophage with replicating phage
  DNA and should not be read as integrated copy number at 37 °C.
- The double-lysogen populations are deliberately not resolved into array
  architecture; read lengths cannot link a marker cassette to one element of a
  tandem array.
- Recurrent insertion sites cannot be separated into insertion hotspots versus
  standing variation with these data.
- Assumption and residual diagnostics: not assessed.

## 8. Files

| File | Contents |
|---|---|
| `timeseries_by_condition.csv` | 11 groups × 4 days, every metric, with n_flasks and min/max |
| `timeseries_per_sample.csv` | the same per library |
| `genotype_concordance.csv` | expected vs observed genotype, all 65 |
| `lesions_per_sample.csv` | five lesion classes per gene per library |
| `lesion_catalogue.csv` | the four point lesions that passed every filter |
| `regression_day8.csv` | 121 day-8 metrics, old vs re-run |
| `read_stats.csv`, `sample_sheet.csv` | library QC and the run definition |
| `fig7`–`fig9` | copy number, marker/lysis, receptor/activator |
