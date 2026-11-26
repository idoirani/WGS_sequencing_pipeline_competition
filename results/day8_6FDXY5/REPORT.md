# Competition of lysis-competent and lysis-deficient λ cI857 lysogens

**Nanopore whole-genome analysis of eleven evolved populations, run 6FDXY5**

---

## 1. Executive summary

Eleven *E. coli* populations carrying λ cI857 prophages were competed for eight days, at
30 °C (repressor active, prophage quiescent) or 37 °C (repressor inactivated, prophage
induced), in defined starting mixtures of a lysis-competent and a lysis-deficient
lysogen. Each population was sequenced to 39–355× with Oxford Nanopore long reads.

Five findings, in order of strength:

1. **Neither competitor won; all four mixtures drifted toward parity.** Final
   lysis-deficient fractions span 0.505–0.780, with both 50/50 flasks ending slightly above
   parity (0.590, 0.571) and both 20/80 flasks moving down toward it (0.780, 0.505). The
   movement is 3–4× larger at 37 °C than at 30 °C, matching where induction occurred. The
   shape is that of frequency-dependent fitness, but with one flask per condition it cannot
   be tested against replicate variance.

2. **Induction at 37 °C roughly doubles prophage DNA per chromosome**, from 1.02–1.51
   copies at 30 °C to 1.93–2.87 at 37 °C — direct confirmation that the temperature
   shift did what it was meant to do.

3. **Receptor loss is the dominant evolutionary response, and it is near-fixation.** At
   37 °C, 72.6–84.3% of reads spanning *lamB* carry a disrupting lesion; at 30 °C the
   figure is 1.0–9.5% and consists entirely of single-base indels in a homopolymer —
   the Nanopore error floor, measured here at 4.9% (19/392 reads). This is the clearest
   signal in the dataset and it was invisible to the assembly-based analysis.

4. **Receptor loss is polyphyletic.** Fifteen distinct *lamB* insertion junctions were
   recovered (11 IS*1*, 3 IS*2*, 1 IS*5*) at fourteen different positions, plus a
   nonsense mutation (Trp99→stop, H7, 19–20%) and a frameshift (codon 333, H11, 73%).
   Independent lesions, not one clone sweeping.

5. **The resident prophage recombined with a cryptic prophage during the experiment.**
   16.4% of H7 reads and 1.2% of H9 reads across the λ tail-fibre locus carry *tfa*
   sequence donated by the Rac cryptic prophage *tfaR* region — the canonical lambdoid
   mosaic hotspot, caught in the act inside a single 8-day experiment. This is the one
   result that speaks directly to the project's recombination hypothesis.

A genome-wide joint variant scan over all 4.7 Mb × 11 samples returned **167 raw sites**
and, after subtracting ancestral differences, requiring sample-specificity, and masking
multi-copy context, **exactly two evolved point mutations — both in *lamB***. Nothing
else in the genome rose to detectable frequency. The scan is therefore an independent
confirmation of finding 3 and a well-controlled negative result everywhere else.

---

## 2. Design and material

Three ancestral lysogens, all sharing the same λ cI857 backbone and differing only at the
right end of the prophage (Figure 1):

| Strain | Prophage | Marker | Lysis |
|---|---|---|---|
| TG5HVC_1 | λ cI857 *bor*::kanR | *aph(3′)-Ia* | competent |
| DD7VVC_2 | λ cI857 *bor*::cmR | *catA1* | competent |
| DD7VVC_1 | λ cI857 Δ(*S*–*R*–*Rz*)::kanR | *aph(3′)-IIa* | **deficient** |

Both lysis-competent markers replace the same 2.06 kb *bor*–*p*78–*p*79 segment; the
lysis-deficient marker replaces the *S*–*R*–*Rz* lysis cassette. The two kanamycin
cassettes are different genes — *aph(3′)-Ia* vs *aph(3′)-IIa* — which is what makes the
three competitors separable by sequence despite two of them sharing a resistance
phenotype. TG5HVC_2 (λ WT, unmarked) served as a fourth reference assembly.

Eleven populations: three pure single-strain controls, two lysis⁺/lysis⁺ mixtures, four
lysis⁺/lysis⁻ mixtures, and two double lysogens, split across the two temperatures.

---

## 3. Scope — what this design can and cannot answer

Two limits are structural, not analytical, and are stated here so no downstream reader
over-reads the tables.

**The double lysogens (H6, H11) cannot be resolved into array architecture.** Resolving
whether two prophages sit in tandem at one *att* site requires a single read that spans a
marker cassette and an internal *att* junction, 20–30 kb apart (Figure 2). At the
realised read-length distributions the expected number of such reads is ≈1 at 30 kb and
falls below 1 beyond that; the internal junction itself (2 kb) is covered dozens of times
over. The libraries detect that a tandem exists but cannot phase a marker onto one. No
copy-number or marker-ratio statement in this report should be read as an architecture
claim.

**A consensus assembly of a mixed population reports the majority allele as though it
were the only allele.** This is the methodological spine of the whole analysis. The
assembly-based pass found *lamB* disrupted in two populations; the read-based pass found
it disrupted in five, at 73–84%. Both passes are shown in the figure recap (§7) precisely
because the discrepancy is instructive. Every frequency in this report is computed by
assigning individual reads to alleles, never by comparing assemblies.

---

## 4. Procedures

### 4.1 Sequencing and quality control

Oxford Nanopore reads for the eleven evolved populations and four ancestral assemblies.
Per-sample yield, length distribution and quality were summarised with a single pass over
each FASTQ (`01_qc.py` → `read_stats.csv`). Reads were mapped with **minimap2** in
`map-ont` mode, sorted and indexed with **samtools**; the exact command was recovered from
the `@PG` line of the existing BAMs so the shipped pipeline reproduces the alignments that
produced these results rather than an approximation of them.

### 4.2 Reference construction

All eleven samples were mapped to a single common reference — the TG5HVC_1 ancestral
assembly — so that every sample sits in one coordinate frame. The 1,900 bp marker cassette
interval was hard-masked (`02_build_ref.py` → `ref_probe.fa`); without the mask, reads
carrying a different cassette soft-clip at the boundary and the locus becomes unmappable
for two of the three competitors. Masked and unmasked references were verified to share
identical contig coordinates so that positions are interchangeable between stages.

### 4.3 Allele assignment at the marker locus

Each read overlapping the cassette window was assigned to exactly one of the three marker
alleles by best alignment identity, with a MAPQ ≥ 20 floor and a requirement of ≥ 200 bp of
flanking anchor sequence (`04_locus_assign.py`). Fractions carry Wilson binomial
confidence intervals. **The internal control is that the three fractions must sum to 1**:
observed sums run 0.963–1.039, and exactly 1.000 in all three pure controls. An earlier
estimator based on sequencing-depth ratios failed this test (sums 0.78–1.20) and was
discarded.

### 4.4 Lysis genotype

Independently of the marker, lysis status was read from the *S*–*R*–*Rz* locus itself:
reads were required to span the junction between the lysis-cassette interval and the
shared right flank, and classified as lysis⁺ (cassette present) or lysis⁻ (deleted)
(`05_lysis_junction.py`). This gives a marker-independent estimate of the same quantity,
and the two agree.

### 4.5 Prophage copy number

Prophage copies per chromosome were estimated as the depth ratio of the 45.6 kb prophage
interval to the chromosomal background. The null was calibrated empirically by applying
the identical estimator to size-matched non-prophage windows drawn from the same BAM,
yielding a per-sample standard deviation (0.088–0.207) against which the 37 °C elevation
is judged (`06_copy_number.py`).

### 4.6 Structural variants and receptor lesions

A read-level scan (`07_sv_scan.py`) collected insertions ≥ 200 bp from CIGAR strings and
clipped-read clusters genome-wide, requiring ≥ 2 independent junction reads. Inserted
sequence was extracted and identified against the ancestral assemblies by **blastn**
(`08_extract_inserts.py`), assigning each to an IS family. In parallel, every read
spanning *lamB* was assigned to exactly one allele class — IS insertion, nonsense,
frameshift, or intact — in a single unified pass (`09_receptor_alleles.py`), so that the
mechanism classes are mutually exclusive and sum to the total disruption frequency.
Insertion junctions were tabulated per position and IS family (`10_is_junctions.py`).

### 4.7 Tail-fibre recombination

Reads spanning the λ *tfa* locus were classified as parental or recombinant by whether the
central segment matched λ *tfa* or the Rac prophage *tfaR* donor, requiring λ sequence on
both flanks so that a recombinant call cannot be produced by a read mapping to Rac itself
(`12_tfa_recomb.py`).

### 4.8 Genome-wide variant scan

All eleven BAMs were called jointly with **bcftools mpileup/call** against the mapping
reference, with a per-file depth cap of 150 and a QUAL ≥ 30 floor
(`14_variant_scan.sh`). The raw call set was then reduced by three filters in order
(`15_variant_filter.py`):

1. **Parental subtraction** — sites at which any ancestor already differs from the
   reference (`13_parental_catalogue.sh`).
2. **Differential requirement** — a site must exceed 0.15 in one sample and fall below
   0.03 in another. A variant present in every sample, including the pure controls, did
   not arise during the experiment.
3. **Repeat masking** — the ±450 bp context must be single-copy in the reference.

Filter 3 is load-bearing and was added only after it was needed. A dense block of 56
linked "variants" in a 660 bp window at ~3,351 kb appeared at 0.22–0.48 in exactly the
five samples containing DD7VVC_1 — a pattern that mimics a lineage-specific sweep.
Probing the window against all four ancestral assemblies showed it to be a 3–4 copy
paralogous family at ~97% identity spanning *tfaD* and an adjacent DNA-packaging gene,
of which DD7VVC_1 alone lacks the 100%-identical copy. Reads from diverged copies
therefore collapse onto the reference copy in precisely those samples. The block is an
artefact of repeat structure, not evolution. The same probe-against-ancestors check had
earlier reclassified an apparent *rssB* indel as a parental difference.

Surviving sites were classified by sample-specificity, which **outranks** homopolymer
context: the Nanopore indel floor is a low fraction spread across many samples, whereas a
real indel is high in one sample and absent from the rest. H11's genuine *lamB* frameshift
sits inside a 5×T run and would be discarded by a naive homopolymer filter.

---

## 5. Results

### 5.1 Sequencing

74,097–606,393 reads per sample (0.18–1.58 Gb), read N50 3.2–5.4 kb, median chromosomal
depth 38.9–354.8×, ≥ 99.99% of reads mapped, Q20 94.0% throughout. H5 (38.9×) and H3
(79.2×) are the depth-limited samples and carry correspondingly wide intervals.

### 5.2 Competition outcome

Change in the lysis-deficient fraction from inoculum to day 8:

| Sample | Condition | T (°C) | lysis⁻ start | lysis⁻ day 8 | Δ |
|---|---|---|---|---|---|
| H4 | 50/50 ⁺/⁻ | 30 | 0.500 | 0.590 | **+0.090** |
| H5 | 20/80 ⁺/⁻ | 30 | 0.800 | 0.780 | −0.020 |
| H8 | 20/80 ⁺/⁻ | 37 | 0.800 | 0.505 | **−0.295** |
| H9 | 50/50 ⁺/⁻ | 37 | 0.500 | 0.571 | +0.071 |

Neither competitor won. The four final fractions span 0.505–0.780 and all four moved
**toward parity**: both flasks started at 50/50 ended slightly above it (0.590, 0.571),
and both started at 20/80 ended at or below their start (0.780, 0.505). The movement is
much larger at 37 °C than at 30 °C (|Δ| 0.295 and 0.071 versus 0.090 and 0.020), which is
the direction expected if induction is what drives the ratio — prophage copy number rises
1.9–2.9× at 37 °C versus 1.0–1.5× at 30 °C (§4.5).

This is the pattern a frequency-dependent fitness function produces, with an equilibrium
somewhere near 0.55 lysis⁻. It is not, however, a tested result: there is one flask per
condition, so there is no replicate variance against which to ask whether convergence on
0.5–0.6 is anything more than four draws from a wide distribution. H8's inoculum is also
the *corrected* condition assignment (§8), so the largest single movement is the one most
contingent on the relabelling being right.

Two internal controls support the frequencies themselves. The pure controls returned
1.000/0.000 exactly and the lysis⁺-only mixtures returned no lysis⁻ reads. And the lysis
genotype measured at the *S*–*R*–*Rz* junction agrees with the lysis status predicted
independently from marker composition to within 0.068 across all eleven samples — two
assays on different loci, in agreement.

### 5.3 Prophage copy number

30 °C: 1.018–1.507 copies per chromosome (null SD 0.088–0.178).
37 °C: 1.932–2.867, except H11 at 0.761.

The 37 °C elevation is 5–10 null SDs and is the expected consequence of cI857
inactivation. H11 is the exception in both directions — lowest copy number and highest
*lamB* disruption — consistent with a population that has largely stopped being a
functional lysogen.

### 5.4 Receptor loss

| Sample | T | *lamB* disrupted (95% CI) | IS | nonsense | frameshift |
|---|---|---|---|---|---|
| H1–H6 | 30 | 0.010–0.095 | 0 | 0 | all |
| H7 | 37 | 0.751 (0.687–0.806) | 0.508 | 0.203 | 0.041 |
| H8 | 37 | 0.726 (0.604–0.821) | 0.710 | 0 | 0.016 |
| H9 | 37 | 0.817 (0.754–0.866) | 0.800 | 0 | 0.017 |
| H10 | 37 | 0.824 (0.698–0.904) | 0.804 | 0 | 0.020 |
| H11 | 37 | 0.843 (0.794–0.882) | 0.109 | 0 | 0.734 |

Every 30 °C sample's apparent disruption is frameshift-class only, at or below the 4.9%
error floor; no IS insertion and no nonsense allele occurs at 30 °C at all. The 37 °C
samples are an order of magnitude above that floor by a mechanism the floor cannot
produce.

Fifteen *lamB* junctions were recovered at fourteen distinct positions spanning codons
18–338 (11 IS*1*, 3 IS*2*, 1 IS*5*). H7 alone carries four independent IS*1* insertions
at codons 45, 72, 98 and 338.

### 5.5 Maltose activator

*malT* was disrupted by IS*5* in four 37 °C samples (H7 0.044, H8 0.133, H9 0.071,
H10 0.160) — but at the **identical base pair**, position 330,600, in all four. One
ancestral or early event rising in parallel, not four independent ones. Loss of *malT*
removes *lamB* expression and so is a second route to the same resistance phenotype.

### 5.6 Tail-fibre recombination

H7 0.164 (0.129–0.206), H9 0.012 (0.005–0.028), zero in the other nine samples. Reads
carry λ *tfa* on both flanks with Rac *tfaR* sequence between, so the call requires a
genuine crossover and cannot arise from mis-mapping to Rac.

### 5.7 Genome-wide scan

167 raw sites → 139 after parental subtraction → 64 differential → 8 single-copy → **2
evolved**:

| Position | Change | Gene | Consequence | Frequency |
|---|---|---|---|---|
| 4,325,749 | GTTTTT→GTTTT | *lamB* | frameshift, codon 333 | H11 0.639, others 0.000 |
| 4,325,797 | C→T | *lamB* | Trp99→stop | H7 0.191, others 0.000 |

Both were already known from the targeted read-level assay (H11 frameshift 0.734, H7
nonsense 0.203) and are here recovered independently, by a different method, at
comparable frequency. The remaining six single-copy sites are four homopolymer
error-floor indels and two sites non-zero in nearly every sample including the pure
controls.

---

## 6. Interpretation

The 37 °C populations were under intense λ predation of their own making: cI857
inactivation induced the resident prophage, free phage accumulated, and the only cells
that survived were those that stopped presenting the receptor. Three independent routes
to that phenotype appear simultaneously — IS insertion into *lamB* (dominant), point
lesions in *lamB*, and IS*5* disruption of the *malT* activator — and together they reach
73–84% of the population in every 37 °C sample. The genome-wide scan found no other locus
under comparable selection, which makes receptor loss not merely the strongest signal but
effectively the only one.

Against that backdrop, the lysis competition settled toward coexistence rather than
exclusion, and it did so from both directions. One mechanistic reading is that receptor
loss, available to both competitors equally and reaching 73–84% in every induced flask,
removes much of the pressure that would otherwise distinguish a lysis-deficient prophage
from a lysis-competent one: once most of the population cannot be infected, whether a
prophage can complete lysis matters less to the fate of the cell carrying it. A second is
genuine frequency dependence — the rarer competitor gains, whichever it is. The data are
consistent with both and establish neither; separating them needs replicate flasks and at
least one more time point.

The tail-fibre recombinants are the most direct observation in the dataset bearing on the
project's central hypothesis. A resident prophage exchanged tail-fibre sequence with a
cryptic prophage elsewhere on the same chromosome, at 16% of the population in one flask,
within a single eight-day experiment. That is prophage mosaicism observed on an
experimental rather than a phylogenetic timescale.

---

## 7. Figure recap

Figures are numbered in the order they were produced. Two are explicitly superseded; they
are retained because the correction they received is itself a result.

### Figure 1 — Ancestral constructs
![Figure 1]({{artifact:art_2a03c336-f9a6-4dcf-b83f-7640922fca0d}})

Map of the three ancestral prophages. **a**, the full 48.5 kb integrated prophage from
*attL* to *attR* with module annotation; all three share the backbone and differ only near
*cos*. **b**, 17–21 kb detail: both lysis⁺ markers replace the same 2.06 kb
*bor*–*p*78–*p*79 segment, while the lysis⁻ marker replaces *S*–*R*–*Rz*. Establishes that
the competitors are distinguishable by sequence and that the two kanamycin cassettes are
different genes.

### Figure 2 — Read geometry
![Figure 2]({{artifact:art_6a10baae-461e-477c-abb0-fc8740b4489e}})

**a**, what a single read must span to resolve a tandem array — the internal *att*
junction is the only sequence unique to a tandem, and linking a marker to it requires
18.7–28.2 kb of contiguous sequence. **b**, expected reads spanning a feature as a
function of required length, per library. The 2 kb junction is covered 60–130×; the
20–30 kb marker-to-junction span falls to ≈1 read. This figure defines the scope limit in
§3 and is the reason no array-architecture claim is made.

### Figure 3 — Assembly-based markers and *lamB* — **superseded**
![Figure 3]({{artifact:art_779e8dac-2707-4480-8fab-28b11804ef6a}})

The assembly-based pass, before FASTQ were available. **a**, *lamB* disrupted in two
populations only (H10 IS*2*, H11 1 bp deletion). **b**, marker copies per assembly —
every population appears dominated by a single competitor lineage. Both panels are
artefacts of consensus assembly: the read-based analysis finds *lamB* disrupted in five
populations at 73–84% (Figure 6) and marker compositions that are genuinely mixed
(Figure 5). Retained as the clearest demonstration of why mixed populations must be
analysed at read level.

### Figure 4 — Depth-ratio marker frequency — **superseded**
![Figure 4]({{artifact:art_fa836074-7b4a-4293-ba2e-50fb330a6773}})

First read-based pass, using sequencing-depth ratios. **a**, marker composition — the
per-sample fractions visibly fail to sum to 1 (0.78–1.20), which is what retired this
estimator in favour of direct read assignment. The two gold arrowheads mark H8 and H10,
whose genotype contradicted their assigned condition; the label transposition was later
confirmed. **b**, sequencing vs day-8 plating ratios, showing ~3-fold compression.
**c**, first prophage copy-number estimate, superseded by the null-calibrated version in
Figure 5.

### Figure 5 — Competition outcome
![Figure 5]({{artifact:art_64546dc1-05ae-47ed-9c68-467d78198491}})

The corrected composition analysis: marker fractions by direct read assignment (summing to
0.963–1.039), lysis genotype from the *S*–*R*–*Rz* junction, observed versus predicted
lysis-deficient fraction, and prophage copy number against its empirically calibrated
null. The quantitative basis for §5.2 and §5.3.

### Figure 6 — Receptor loss
![Figure 6]({{artifact:art_a3fdcb4b-594b-46e1-825e-26c775a60b13}})

**a**, *lamB* disruption frequency per sample with Wilson intervals, 30 °C versus 37 °C,
with the measured error floor marked. **b**, mechanism composition — IS, nonsense,
frameshift, intact — showing H7 and H11 reaching the same endpoint by different routes.
**c**, *malT* disruption. **d**, *tfa* recombinant frequency. The quantitative basis for
§5.4–5.6.

---

## 8. Corrections and open questions

**Applied correction.** H8 and H10 carry transposed condition labels. The genotype
contradiction was flagged from the marker composition (Figure 4a) and subsequently
confirmed. Every table carries the corrected condition with the original label retained in
brackets, so the correction remains auditable.

**Unsettled — hotspot versus standing variation.** The four *malT* IS*5* insertions at the
identical base pair, and the three IS*1* insertions shared at codon 338 across H7, H8 and
H9, admit two readings: an insertion hotspot hit repeatedly and independently, or a single
ancestral insertion present at low frequency in the inoculum and enriched in parallel.
The data cannot distinguish these. Deep sequencing of the ancestral stocks would settle
it directly.

**Not assessed.** Statistical significance of the competition deviations against
replicate-level variance; there is one flask per condition. Free-phage titres. Whether
receptor-loss mutants retain the prophage.

---

## 9. Reproducing this analysis

The complete pipeline is in `pipeline/`, with `config.py` as the single edit point for
paths, sample names, coordinates and thresholds, and `run_all.sh` executing the stages in
order. `pipeline/README.md` documents each stage and its outputs; `TOOLS.md` records every
tool and version that produced these results. Four stages were re-run after refactoring
and reproduce byte-identical output.

To run on new data: place FASTQ in `fastq/<sample>.fastq.gz`, edit `SAMPLES`, `TEMP` and
the condition dictionaries in `config.py`, and run `bash pipeline/run_all.sh`.

---

## Appendix — earlier comparative-genomics work in this project

This report covers the competition experiment only. The project also contains a separate
comparative-genomics analysis of prophage recombination and selection, whose figures are
indexed here for completeness: prophage inventories and defectiveness profiles
(`fig_prophage_inventory.png`, `fig_mg1655wt_prophages.png`), panel QC and phylogroup
annotation (`fig1_genome_panel_qc.png`, `fig2_phylogroups_annotation.png`), the prophage
landscape and clustering appendix (`fig3_prophage_landscape.png`,
`fig_clustering_appendix.png`), linkage-disequilibrium decay and its distance-axis
diagnosis (`fig4_recombination.png`, `fig4_recombination_reconciled.png`,
`fig_ld_axis_diagnosis.png`, `fig_axis_schematic.png`, `fig_recombination_tests.png`),
dN/dS by gene group and estimator reconciliation (`fig5_selection.png`,
`fig5_selection_reconciled.png`, `fig_dnds_by_group.png`, `fig_dnds_reconciliation.png`,
`fig6_estimator_comparison.png`), the P22 comparison against the original draft
(`fig_p22_vs_document.png`) and the combined summary (`fig_summary_both_panels.png`).
That analysis has its own `SUMMARY.md`.
