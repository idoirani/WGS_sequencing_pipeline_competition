# Tools

Every version below is the one that produced the results in `REPORT.md`, read
from the installed binaries rather than from a specification file.

## Environment

All tools are **native arm64** on Apple Silicon. This differs from the earlier
comparative-genomics work in this project (geNomad / Pharokka / Roary / Mash),
which had to run under `osx-64` Rosetta 2 emulation. Nothing in the competition
pipeline needs emulation, and nothing in it needs a GPU.

```bash
conda create -y -n phage-evo -c conda-forge -c bioconda \
  python=3.11 minimap2 samtools bcftools seqkit mummer4 blast biopython \
  numpy scipy matplotlib
conda activate phage-evo
```

or equivalently `bash pipeline/00_setup_env.sh phage-evo`.

## Versions

| Tool | Version | Used for |
|---|---|---|
| Python | 3.11.16 | all analysis scripts |
| minimap2 | 2.31-r1302 | ONT read alignment (`-ax map-ont --secondary=no`) |
| samtools | 1.24 | sorting, indexing, region extraction, depth |
| bcftools | 1.24 | joint pileup, variant calling, filtering |
| seqkit | 2.13.0 | read-length and quality statistics |
| MUMmer (nucmer) | 4.0.1 | parental assembly comparison |
| BLAST+ (blastn) | 2.17.0+ | identifying inserted elements, uniqueness checks |
| Biopython | 1.88 | FASTA/GenBank parsing, reference construction |
| NumPy | 2.4.6 | array arithmetic, bootstrap resampling |
| SciPy | 1.17.1 | binomial confidence intervals |
| Matplotlib | 3.11.1 | figures |

## What each tool is load-bearing for

**minimap2** — `--secondary=no` is not cosmetic. Every assay in this pipeline
counts primary alignments only; if secondary alignments were kept, a single read
could be counted against two marker alleles at once and the allele fractions
would not sum to 1.

**samtools** — `-F 0x904 -q 20` (drop unmapped, secondary and supplementary;
MAPQ ≥ 20) is the read filter used by every counting step. Changing it changes
every reported fraction.

**bcftools** — the per-file depth cap (`-d`) is the one parameter that decides
whether the genome-wide scan finishes at all. See `pipeline/14_variant_scan.sh`.

**nucmer** — used only to compare the four parental assemblies to each other, so
that pre-existing differences between the competing strains are subtracted from
the evolved-variant call set rather than being reported as mutations.

**blastn** — used to identify extracted insertion sequences against the strain's
own chromosome, and to confirm that every interval used as a denominator is
single-copy before it is used as one. The latter use is load-bearing: the
cryptic-prophage *tfa*/*tfaD*/packaging genes form a 3–4 copy family at ~97%
identity, and reads from a diverged copy collapsing onto the reference copy
generate dense blocks of convincing-looking, strictly lineage-correlated false
variants. `pipeline/15_variant_filter.py` blasts the context of every candidate
site and discards anything with more than one hit.

**seqkit** — also used in `13_parental_catalogue.sh` to pull a single contig out
of a reference FASTA. `samtools faidx` is the obvious tool for that and cannot be
used here: it insists on writing a `.fai` beside the input, and the ancestral
assemblies live on a read-only mount.

## Not used, and why

- **No assembler.** These are mixed populations. A consensus assembly of a
  mixed population reports the majority allele as if it were the only allele,
  which is exactly the quantity this experiment is measuring. All calls are made
  from reads against a fixed reference.
- **No variant caller tuned for diploids/somatic data.** Allele fractions are
  counted directly from reads at defined loci, where the denominator can be
  audited, rather than inferred from a caller's genotype likelihoods.
- **No phasing or strain deconvolution.** The double-lysogen populations are
  deliberately not resolved into constituent genomes.
