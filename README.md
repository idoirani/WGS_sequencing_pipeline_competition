# WGS sequencing pipeline — lysogen competition

Quantifies competing lysogen strains **inside a mixed bacterial population**
from Oxford Nanopore whole-genome reads: marker-cassette allele frequencies,
lysis genotype, prophage copy number, phage-receptor (*lamB*) and activator
(*malT*) disruption, tail-fibre recombinants, and genome-wide evolved variants.

Built for λ *cI857* lysogens competing lysis-competent against lysis-deficient
in serial passage, but the structure holds for any experiment where two or more
engineered strains are mixed and you need their ratio over time.

---

## The rule that shapes the whole design

**Never assemble a mixed population.** A consensus assembly reports the majority
allele as if it were the only allele — which is exactly the quantity being
measured. In the day-8 dataset this pipeline was built on, assemblies called one
cassette per sample and hid the minority competitor in *all eight* mixtures, and
reported the receptor disrupted in 2 of 11 populations where the reads show it
disrupted in 5 at 73–84 %.

Every number here is counted from reads against a fixed reference, with a
denominator you can audit.

---

## Quick start

```bash
bash pipeline/00_setup_env.sh phage-evo     # conda env from environment.yml
conda activate phage-evo

$EDITOR pipeline/config_inputs.py           # where this experiment's data is
$EDITOR pipeline/config.py                  # the biological system (rarely changes)

bash pipeline/run_all.sh                    # all stages
bash pipeline/run_all.sh 04                 # resume from stage 04
KEEP_BAM=1 bash pipeline/run_all.sh         # keep whole-genome BAMs (stage 14)
```

Mapping is resumable: a library whose region BAM and depth table already exist
is skipped, so an interrupted run costs nothing to restart.

## The two config files

| File | Holds | Changes |
|---|---|---|
| `pipeline/config_inputs.py` | paths to the plate record and the read archives | every experiment |
| `pipeline/config.py` | reference, coordinates, thresholds | only if the strains change |

`config.py` carries no sample names: the run is defined by `sample_sheet.csv`,
which stage 00 builds from the lab's own plate record.

## Sample sheet

Stage `00_sample_sheet.py` reads the plate map (one block per plate, rows =
passage day, columns = condition) and joins it to the read archives. A library
named `v1B3` is plate 1, row B (day 2), column 3 — so replicate, day and
condition follow from the library name and nothing is typed in by hand.

Columns consumed by the pipeline:

| Column | Meaning |
|---|---|
| `sample` | library name; every output is keyed on it |
| `fastq` | path, glob, or `zip:<archive>!<member>` to stream from a zip |
| `temp_C`, `day`, `rep` | experiment coordinates |
| `condition`, `group` | the lab's own label, and a key shared across days |
| `control_allele` | non-empty for pure single-strain controls — the allele that library must show at 1.0 |
| `lysis_plus_only` | 1 where the inoculum was entirely lysis-competent |
| `use` | 0 excludes a failed library without deleting the row |

Reads are streamed from the archive member-by-member; a 47 GB zip is never
unpacked.

**Controls are not optional.** Near-identical alleles cross-assign, and the
pure controls are the rows of the confusion matrix that stage 05 inverts. With
no control for an allele, stage 05 stops rather than reporting biased fractions.

## Stages

| Stage | Produces |
|---|---|
| `00_setup_env.sh` | conda environment |
| `00_sample_sheet.py` | `sample_sheet.csv` from the plate record + archives |
| `02_build_ref.py` | masked backbone + one decoy contig per marker allele |
| `02b_build_locus_ref.py` | four-variant locus reference, byte-identical anchors |
| `03_map.sh` | region BAM + 1 kb depth table + flagstat per library |
| `01_qc.py` | `read_stats.csv` |
| `04_locus_assign.py` | marker allele fractions by direct read assignment |
| `05_lysis_junction.py` | lysis genotype from the *S–R–Rz* junction; unmixing |
| `06_copy_number.py` | prophage copies/chromosome, faux-prophage null |
| `07_sv_scan.py` | structural events vs the matching ancestor |
| `08_extract_inserts.py` | inserted-element sequences pulled from reads |
| `09_receptor_alleles.py` | *lamB* lesion classes, one class per read |
| `10_is_junctions.py` | IS family and exact junction per insertion |
| `11_activator.py` | *malT* lesions |
| `12_tfa_recomb.py` | tail-fibre recombinant frequency |
| `13_parental_catalogue.sh` | sites that differ between the ancestors |
| `14_variant_scan.sh` | genome-wide joint pileup and call |
| `15_variant_filter.py` | parental / differential / **repeat** filters, then classify |
| `30_tables.py` | result CSVs |

## Five things that cost hours if skipped

1. **Map with `--secondary=no`; count with `-F 0x904 -q 20`.** Otherwise one
   read is credited to two alleles and the fractions stop summing to 1.
2. **Do not use depth ratios at an engineered locus.** A coverage shadow around
   the mask plus an induction-scaled gradient across the prophage left allele
   budgets at 0.77–1.04. Direct read assignment against identical flanks closed
   them to 0.963–1.039. Stage 04 is the method of record.
3. **Unmix against pure controls** (13.2 % cross-assignment measured here).
4. **Mask repeats before believing any variant.** The cryptic-prophage
   *tfa*/packaging family is 3–4 copies at ~97 % identity; reads from a diverged
   copy collapse onto the reference copy and manufacture a dense block of
   linked, strictly lineage-correlated false variants — a convincing fake sweep.
   Stage 15 blasts ±450 bp of context and keeps only single-hit sites.
5. **Sample-specificity outranks homopolymer context** when classifying indels.
   The ONT indel floor is a low fraction spread over many samples; a real indel
   is high in one and absent from the rest. Reversing that order discards
   genuine frameshifts sitting in mononucleotide runs.

## Failure modes with non-obvious causes

- **Variant scan returns nothing, instantly or after hours.** The BAMs are
  mapped to the *masked* reference from stage 02, not the plain ancestor FASTA.
  Point `bcftools mpileup` at `MAP_FASTA`. Check `MPILEUP_D` too — an unbounded
  depth cap stalls a multi-sample pileup indefinitely.
- **`samtools faidx` fails on the reference.** It writes its `.fai` beside the
  input, and reference assemblies usually sit on a read-only mount. Use
  `seqkit grep` (stage 13 does).
- **Region BAM will not index.** Passing several regions to `samtools view`
  emits reads in *region* order, not coordinate order, and emits a read twice
  where two regions overlap. Stage 03 uses a merged, sorted BED with `-L`.

## Reporting conventions

Every number entering a figure caption or a report is computed from a saved
table and read back before it stands. Distinct lineages and samples stay
individually identifiable in figures unless the caption states a reason to pool
them. Where a label transposition is confirmed, the correction is applied with
the original label retained so it stays auditable. With one flask per condition,
movement is reported as a pattern with the absence of replicate variance stated
explicitly.

## Requirements

`environment.yml` pins the versions that produced the published results
(minimap2 2.31, samtools/bcftools 1.24, seqkit 2.13.0, MUMmer 4.0.1, BLAST+
2.17.0, Python 3.11). `requirements.txt` covers the Python layer only — the
external binaries must come from conda. Native arm64 on Apple Silicon; no
Rosetta anywhere in this pipeline.
