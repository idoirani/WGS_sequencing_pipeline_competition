# Analysis record — CI857 Δlysis vs lysis-competent lysogen competition

**Experiment:** `CI857 deltalysis vs lysis competant 10212026`
**Runs analysed:** day 8 only (run 6FDXY5, 11 libraries) → full time course (run KM89QV + 6FDXY5, 65 libraries, days 1–8)
**Analysis dates:** 11–12 Nov 2026, 25 Nov 2026

This file records what was asked, what was done, and what had to be corrected.
The science is in `results/*/REPORT*.md`; the tools are in `TOOLS.md`; the
pipeline and how to re-run it are in `README.md`.

---

## Part 1 — day 8 (11 libraries)

### Asked
1. Resolve the sample ↔ condition mapping, including libraries whose genotype
   contradicts their assigned condition.
2. Pointer to the FASTQ location.
3. *"what about the other mutations? i.e., did any of the strains lose lamB or
   other interesting mutations?"*
4. Executive summary, figure recap, tool markdown, paper-level procedures,
   re-runnable from scratch.
5. Package the pipeline as a reusable skill.
6. Save figures and a conversation record to disk.

### Found
- The four ancestors are **isogenic** on the chromosome (100.000 % over 4.68 Mb);
  the engineered cassette is the only lineage marker. The two kanR competitors
  carry **different** phosphotransferases — aph(3')-Ia (lysis⁺) vs aph(3')-IIa
  (Δlysis), ~65 % identical — so they are sequence-distinguishable even though
  plating cannot separate them.
- **The central methodological result:** a consensus assembly of a mixed
  population reports the majority allele as if it were the only allele.
  Assemblies hid the minority competitor in all eight mixtures and called the
  receptor disrupted in 2 of 11 populations where reads show 5 at 73–84 %.
  Everything was recomputed from reads against a fixed reference.
- *lamB* disrupted in every 37 °C sample: 15 distinct IS junctions (11 IS1,
  3 IS2, 1 IS5) plus a Trp99→stop and a codon-115 frameshift. A *malT* IS5
  insertion at the identical base pair in four samples. Tail-fibre recombinants
  with Rac-prophage donor sequence at 16.4 % in one sample.
- Genome-wide scan: 167 raw sites → **exactly two evolved point mutations in the
  whole chromosome, both in *lamB***, reproducing the targeted assay. A dense
  apparent "sweep" at ~3,351 kb was shown to be a collapsed-paralogue artefact.
- An H8↔H10 label transposition, confirmed by the lab.

---

## Part 2 — full time course (65 libraries, days 1–8)

### Asked
> *"I want to run the pipeline on all of these conditions … find as an output the
> average prophage copy number as a function of time for each condition group,
> the cmR/kanR ratio, and lysis deficient/competant ratio, and the fraction of
> lamB deletions and malT deletions as a function of time per condition group."*

Plus: track the pipeline in a git repository, include a README and dependency
list, and keep the finished outputs with the code.

### Found
1. **Induction produces a burst of phage DNA that collapses.** 37 °C group means
   peak at 4.1–26.3 phage genomes per chromosome at day 1–2 (highest single
   library 49.7) and fall to 0.8–2.9 by day 8. 30 °C never leaves 0.97–2.9.
2. **Receptor loss is late, near-total, and independent of the competition.**
   *lamB* disruption (37 °C means): 0.015 → 0.021 → 0.391 → **0.826** across days
   1, 2, 4, 8. It reaches 0.76–0.87 in **every** 37 °C group at day 8, including
   those with no lysis-deficient competitor. At 30 °C it never exceeds 0.044.
3. **A convergent 25 bp insertion at *malT* codon 341** in 12 of 65 libraries
   above 5 % (max 0.284), only at 37 °C, in both replicates, and exactly 0.000 in
   all 17 pure 30 °C controls. *malT* activates *lamB*, so this is a second,
   earlier route to the same resistance phenotype. *malT* disruption peaks at day
   4 (mean 0.452) and falls to 0.268 by day 8 while *lamB* keeps climbing.
4. **The competition does not resolve.** In ⁺:⁻ mixtures the lysis-deficient
   fraction oscillates around its inoculum with swings as large as the difference
   between the two inocula. The day-8-only "convergence toward parity" reading
   from Part 1 is **not** supported once the intermediate days are included.
5. **Seven libraries disagree with their label**, reported and not relabelled —
   notably a pure control carrying 16.7 % of another cassette, and an apparent
   transposed pair. The lab's own plate record already lists columns 8 and 10
   swapped in the day-8 row, independently confirming the Part 1 transposition.

---

## Corrections made along the way

Recorded because each would have changed a reported number.

| What | Effect if uncorrected |
|---|---|
| `f_lysis_plus_predicted` read as the inoculum (Part 1) | reversed the reported competition outcome |
| *malT* strand `+` in config; the annotation says `−` | every *malT* codon number mirrored about the gene's middle |
| Confusion matrix pooled across plate rows | the IIa→wild-type leak varies 0.055–0.139 between rows; matching tightened allele budgets from 0.944–1.091 to 0.963–1.077 |
| A control that is not pure used as a matrix row | one declared-pure control carries 16.7 % of another cassette, which would bias all 65 libraries |
| Indels keyed by length | ONT reports one true 25 bp insertion as 24–27 bp, splitting the real *malT* lesion into four sub-threshold sites |
| Relative point-lesion error floor | a 1 bp deletion at 5.8 % in every control produced 9–43 % false disruption in pure controls |
| Group key stripping punctuation | pooled `+/-` with `+/+`, and `lysis+` with `lysis-` |
| Region extraction with multiple `samtools view` regions | unindexable BAM, and reads counted twice where regions overlap |

Before any new library was trusted, the repaired pipeline was re-run on the 11
day-8 libraries: **121 of 121 metrics reproduce the saved day-8 tables exactly**
(`results/timeseries_2026-11/tables/regression_day8.csv`).

---

## To re-run

```bash
bash pipeline/00_setup_env.sh phage-evo && conda activate phage-evo
$EDITOR pipeline/config_inputs.py      # where this experiment's data is
bash pipeline/run_all.sh
```

Also available as the skill `lysogen-competition-ont`:
`deploy_competition_pipeline()` writes `pipeline/` into any workspace.

## Known limitations

- No group-day cell has more than two flasks; nothing is tested for significance.
- Day 1 is replicate 2 only and day 8 is replicate 1 only, so day 1 vs day 8 is
  also a between-replicate comparison.
- Copy number after induction mixes integrated prophage with replicating phage
  DNA and is not integrated copy number at 37 °C.
- Double-lysogen populations are deliberately not resolved into array
  architecture; read lengths cannot link a cassette to one element of a tandem array.
- Recurrent insertion sites cannot be separated into hotspots versus standing
  variation with these data.
