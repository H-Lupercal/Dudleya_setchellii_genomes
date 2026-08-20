# Dudleya Organelle Alignment Pipeline

This directory contains the custom cpDNA/mtDNA workflow for the downloaded
Dudleya whole-genome FASTQ data. The workflow is organised into small, auditable
steps so that sample identity and QC are established before bulk alignment.

The authoritative, ordered stage index is [`PROCESS.md`](PROCESS.md). Stages are
identified by their `results/NN_.../` directory number, and the section headings
below use those canonical stage numbers.

## Stage 00: Manifest And Preflight Validation

Stage 00 scans the downloaded FASTQ files and writes a sample manifest. It does
not trim reads, align reads, call variants, or build cpDNA/mtDNA consensus
sequences.

The purpose of this step is to answer these questions first:

- Which FASTQ files were downloaded?
- Which files are `R1` and which files are `R2`?
- Which `R1` and `R2` files belong to the same biological sample?
- Which sequencing batch did each sample come from?
- Which filename convention does each sample use?
- Which samples have population-code metadata available?
- Which samples are safe candidates for the pilot organelle alignment?
- Which samples must be excluded from the primary paired-end analysis?

## Naming Profiles

The parser recognizes three filename profiles that share the Illumina suffix:

```text
_S<sampleNumber>_L<lane>_R<1-or-2>_001.fastq.gz
```

### Main Standard Dataset

Example:

```text
CY_RED_LP_202_Du-561_S192_L005_R1_001.fastq.gz
```

Parsed as:

- `sample_id`: `CY_RED_LP_202_Du-561`
- `popcode`: `CY_RED`
- `lp_id`: `LP_202`
- `du_id`: `Du-561`
- `naming_profile`: `main_standard`

Population metadata is attached from:

```text
genomicsDrive_data_dump/QB3.Berkeley.251217/Dudleya DNAx - Population Codes.csv
```

### Initial DU-Dash Dataset

Example:

```text
DU-4A_S68_L008_R1_001.fastq.gz
```

Parsed as:

- `sample_id`: `DU-4A`
- `du_id`: `DU-4A`
- `naming_profile`: `initial_du_dash`

These samples can be aligned, but population metadata remains unresolved until
a manual lookup table is added.

### Initial DU/LP Dataset

Example:

```text
DU014LP012_S4_L005_R1_001.fastq.gz
```

Parsed as:

- `sample_id`: `DU014LP012`
- `du_id`: `DU014`
- `lp_id`: `LP012`
- `naming_profile`: `initial_du_lp`

These samples can be aligned, but population metadata remains unresolved until
a manual lookup table is added.

## Run Stage 00

From the repository root:

```bash
python3 dudleya_organelle_alignment_pipeline/scripts/build_sample_manifest.py
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/00_manifest/samples.tsv
dudleya_organelle_alignment_pipeline/results/00_manifest/analysis_samples.tsv
dudleya_organelle_alignment_pipeline/results/00_manifest/excluded_samples.tsv
dudleya_organelle_alignment_pipeline/results/00_manifest/pairing_report.tsv
dudleya_organelle_alignment_pipeline/results/00_manifest/preflight_summary.md
```

To reproduce the current preflight state from a clean checkout with downloaded
FASTQs in `genomicsDrive_data_dump/`, run:

```bash
python3 -m unittest dudleya_organelle_alignment_pipeline.tests.test_manifest -v
python3 dudleya_organelle_alignment_pipeline/scripts/build_sample_manifest.py
```

## Output Files

`samples.tsv` has one row per discovered biological sample. It records the
batch, naming profile, sample IDs, population metadata where available, paired
FASTQ paths, pair status, primary-analysis status, and a note explaining the
analysis decision.

`analysis_samples.tsv` is the authoritative input table for the primary
paired-end cpDNA/mtDNA alignment. Future alignment scripts should read this
file, not `samples.tsv`, when running the main analysis.

`excluded_samples.tsv` records samples excluded from the primary paired-end
alignment and the reason for exclusion.

`pairing_report.tsv` records problems that should be reviewed before alignment,
such as missing `R1`, missing `R2`, uneven read counts, or unparsed FASTQ names.

`preflight_summary.md` is the human-readable overview of the run.

## Missing-Mate Policy

Two samples were manually verified as missing their mate FASTQ:

- `ABAB_MAD_LP_225_Du-592`: missing `R1`
- `QUI1_LP_256_Du-655`: missing `R2`

These samples remain in `samples.tsv` for the audit trail, but are excluded
from `analysis_samples.tsv` and therefore from the primary paired-end cpDNA/mtDNA
alignment.

If either sample is ever aligned by itself as a single-end case, that run must
be written up separately as an individual/sensitivity analysis. Do not mix
single-end results from these samples into the primary paired-end alignment,
PCA, tree, Fst, or structure/admixture-style outputs.

## Relationship To The Published Conservation-Genomics Pipeline

A published Dudleya conservation-genomics pipeline (Hackstadt,
https://github.com/evanhackstadt/dudleya) provides the reference pattern for
general whole-genome-sequencing QC: sample-table-driven processing, read QC,
aggregate QC summaries, and organized downstream outputs. This organelle pipeline
starts with a custom manifest because the cpDNA/mtDNA task uses a different
reference structure and needs separate chloroplast and mitochondrial outputs.

After this manifest step, subsequent steps reuse that QC organization while
replacing the biological core with organelle-specific mapping, coverage,
consensus, and cpDNA/mtDNA alignment generation.

## Stage 01: Reference And Pilot Preflight

Stage 01 validates the combined cpDNA/mtDNA reference, records whether required
tools are installed, creates reference indexes only when those tools are
available, and writes a representative pilot sample table.

Run from the repository root:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/prepare_reference_and_pilot.py
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/01_reference_pilot/reference_checks.tsv
dudleya_organelle_alignment_pipeline/results/01_reference_pilot/tool_checks.tsv
dudleya_organelle_alignment_pipeline/results/01_reference_pilot/index_checks.tsv
dudleya_organelle_alignment_pipeline/results/01_reference_pilot/pilot_samples.tsv
dudleya_organelle_alignment_pipeline/results/01_reference_pilot/reference_pilot_summary.md
```

Current result:

```text
Reference records checked: 2
Pilot samples written: 15
```

The combined reference passed the expected record checks:

- `chloroplast`: 150274 bp
- `mitochondria`: 243359 bp

The current pilot table contains 15 complete paired-end samples:

- 2 unresolved initial-batch representatives.
- 5 *D. cymosa* samples.
- 3 *D. abramsii* samples.
- 5 *D. setchellii* samples.

The two manually verified missing-mate samples are not eligible for the pilot
set because `pilot_samples.tsv` is selected from `analysis_samples.tsv`.

## Local Tool Environment

This repo currently has a local, git-ignored tool environment at:

```text
.tools/bioconda-env/
```

Use it by prefixing tool-dependent commands with:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" <command>
```

The project tool path contains the executables needed for the organelle
alignment, variant, tree, PCA, Fst, and admixture/structure-style stages:

```text
bwa 0.7.19-r1273
samtools 1.23.1
bcftools 1.23.1
fastp 1.3.5
FastQC 0.12.1
MultiQC 1.35
IQ-TREE 3.1.2
PLINK 1.9
ADMIXTURE 1.3.0
VCFtools 0.1.17
BEDTools 2.31.1
Rscript 4.5.3
matplotlib 3.11.0
pandas 3.0.3
scikit-learn 1.9.0
Biopython 1.87
seaborn 0.13.2
ggplot2 4.0.3
ape 5.8.1
patchwork 1.3.2
```

The current machine-readable tool audit is recorded in:

```text
dudleya_organelle_alignment_pipeline/results/13_tool_audit/primary.tool_audit.tsv
dudleya_organelle_alignment_pipeline/results/13_tool_audit/primary.tool_audit_report.md
```

If the local environment needs to be recreated, use:

```bash
micromamba create -y \
  -p .tools/bioconda-env \
  -f dudleya_organelle_alignment_pipeline/environment.yml
```

Stage 01 has created the `samtools faidx` and `bwa index` files for:

```text
dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa
```

To reproduce Stages 00 and 01:

```bash
python3 -m unittest \
  dudleya_organelle_alignment_pipeline.tests.test_manifest \
  dudleya_organelle_alignment_pipeline.tests.test_prepare_reference_and_pilot -v
python3 dudleya_organelle_alignment_pipeline/scripts/build_sample_manifest.py
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/prepare_reference_and_pilot.py
```

## Stage 02: Pilot Organelle Alignment

Stage 02 aligns the 15 pilot samples from `pilot_samples.tsv` to the combined
cpDNA/mtDNA reference. It keeps mapped read records, writes sorted/indexed BAMs,
and summarizes cpDNA and mtDNA mapped read counts, input mapping fraction, mean
depth, and breadth at `>=1x`, `>=5x`, and `>=10x`.

It does not call variants, create consensus FASTAs, make final all-sample
alignments, or run PCA/tree/Fst/admixture analyses.

Run from the repository root:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_pilot_alignment.py
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/02_pilot_alignment/pilot_alignment_sample_summary.tsv
dudleya_organelle_alignment_pipeline/results/02_pilot_alignment/pilot_alignment_by_organelle.tsv
dudleya_organelle_alignment_pipeline/results/02_pilot_alignment/pilot_alignment_report.md
dudleya_organelle_alignment_pipeline/results/02_pilot_alignment/commands.tsv
dudleya_organelle_alignment_pipeline/results/02_pilot_alignment/bam/
dudleya_organelle_alignment_pipeline/results/02_pilot_alignment/qc/
dudleya_organelle_alignment_pipeline/results/02_pilot_alignment/logs/
```

The BAM, depth, and log files are generated analysis artifacts and are ignored
by git. The top-level Stage 02 summaries are small enough to keep as the pilot QC
record.

Current pilot result:

```text
Pilot samples summarized: 15
Sample-by-organelle rows: 30
Total cpDNA+mtDNA mapped read records: 46584669
Median input organelle mapping fraction: 0.144205
Median chloroplast breadth >=1x: 0.999993
Median mitochondrial breadth >=1x: 0.960445
```

After correcting the `samtools depth` quality flags, only the tiny
`ABAB_MAD_LP_222_Du-589` pilot sample has an initial QC note. The key mtDNA
issue is now repeat/ambiguity handling: a high-MAPQ-only check shows much lower
unique-placement breadth than permissive MAPQ depth.

The focused mtDNA investigation outputs are:

```text
dudleya_organelle_alignment_pipeline/results/03_mtdna_investigation/mtdna_investigation_report.md
dudleya_organelle_alignment_pipeline/results/03_mtdna_investigation/mtdna_depth_filter_comparison.tsv
dudleya_organelle_alignment_pipeline/results/03_mtdna_investigation/mtdna_high_mapq_consensus_intervals.tsv
dudleya_organelle_alignment_pipeline/results/03_mtdna_investigation/mtdna_high_mapq_coverage_distribution.tsv
```

The focused cpDNA verification outputs are:

```text
dudleya_organelle_alignment_pipeline/results/04_cpdna_investigation/cpdna_verification_report.md
dudleya_organelle_alignment_pipeline/results/04_cpdna_investigation/cpdna_depth_filter_comparison.tsv
dudleya_organelle_alignment_pipeline/results/04_cpdna_investigation/cpdna_self_repeat_intervals.tsv
dudleya_organelle_alignment_pipeline/results/04_cpdna_investigation/cpdna_repeat_overlap_summary.tsv
```

The cpDNA verification supports moving forward with all-sample chloroplast
processing. The main caution is the expected chloroplast inverted repeat:
`82091-107826` and `124539-150274` in the normalized reference.

## Stage 05: Analysis Masks And Tracks

Stage 05 converts the cpDNA and mtDNA verification findings into the exact tracks
that Stage 06 must use. It does not align reads, call variants, create consensus
FASTAs, or run population-genetic analyses.

Run from the repository root:

```bash
python3 dudleya_organelle_alignment_pipeline/scripts/build_analysis_masks.py
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/05_analysis_masks/analysis_tracks.tsv
dudleya_organelle_alignment_pipeline/results/05_analysis_masks/analysis_regions.tsv
dudleya_organelle_alignment_pipeline/results/05_analysis_masks/mask_summary.md
dudleya_organelle_alignment_pipeline/results/05_analysis_masks/cpdna_full_coverage_regions.bed
dudleya_organelle_alignment_pipeline/results/05_analysis_masks/cpdna_ir_regions.bed
dudleya_organelle_alignment_pipeline/results/05_analysis_masks/cpdna_duplicate_ir_copy_mask.bed
dudleya_organelle_alignment_pipeline/results/05_analysis_masks/cpdna_population_sites.bed
dudleya_organelle_alignment_pipeline/results/05_analysis_masks/mtdna_permissive_coverage_regions.bed
dudleya_organelle_alignment_pipeline/results/05_analysis_masks/mtdna_high_confidence_unique_regions.bed
```

Current Stage 05 decisions:

- BED files are 0-based, half-open.
- `analysis_regions.tsv` records the same intervals as 1-based inclusive
  coordinates plus their BED coordinates.
- cpDNA sample-level coverage QC uses the full 150,274 bp chloroplast
  reference.
- cpDNA PCA, Fst, tree, and clustering inputs use
  `cpdna_population_sites.bed`, which retains 124,538 bp by keeping one IR copy
  and excluding the duplicate IR copy at `124539-150274`.
- mtDNA sample-level QC stays in the project as a whole-reference permissive
  coverage track.
- mtDNA variant calling and population-genetic inputs use
  `mtdna_high_confidence_unique_regions.bed`, currently two high-MAPQ consensus
  intervals totaling 44,930 bp.

`analysis_tracks.tsv` is the authoritative machine-readable contract for Stage 06:
coverage-QC tracks are not interchangeable with population-genetic tracks.

## Stage 06: All-Sample Organelle Alignment

Stage 06 maps every primary paired-end sample from `analysis_samples.tsv` to the
combined cpDNA/mtDNA reference and summarizes coverage by organelle and by the
Stage 05 analysis tracks.

It does not call variants, create consensus FASTAs, make final alignments, or
run PCA/tree/Fst/admixture analyses.

Run from the repository root:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_all_sample_alignment.py
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/all_sample_alignment_sample_summary.tsv
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/all_sample_alignment_by_organelle.tsv
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/all_sample_alignment_by_track.tsv
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/all_sample_alignment_report.md
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/downstream_sample_qc_decisions.tsv
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/downstream_sample_qc_decisions.md
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/commands.tsv
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/bam/
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/qc/
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/logs/
```

The BAM, depth, and log files are generated analysis artifacts and are ignored
by git. The top-level Stage 06 summaries are small enough to keep as the
all-sample QC record.

The command is resumable. If an output BAM and its QC files already exist for a
sample, Stage 06 reuses them unless `--force` or `--refresh-qc` is passed.

After reviewing Stage 06 QC, the downstream primary analysis set excludes:

- `ABAB_MAD_LP_222_Du-589`
- `CY_HUN_LP_265_Du-684`
- `CY_RED_LP_202_Du-561`

These samples are ignored downstream because they are the three lowest-input
samples in the run and failed one or both organelle coverage screens. The
machine-readable decision record is
`results/06_all_sample_alignment/downstream_sample_qc_decisions.tsv`.

## Stage 07: Downstream Sample Set

Stage 07 converts the Stage 06 QC decisions and upstream missing-mate exclusions
into the exact sample tables that variant calling and later population-genetic
steps must use.

Run from the repository root:

```bash
python3 dudleya_organelle_alignment_pipeline/scripts/build_downstream_sample_set.py
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/included_samples.tsv
dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/excluded_samples.tsv
dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/downstream_sample_set_report.md
```

The primary downstream sample set should contain 275 samples.

## Stage 08: Haploid Variant Calling

Stage 08 uses the Stage 07 included sample set and calls raw haploid variants
separately for cpDNA and mtDNA. It restricts cpDNA to the IR-aware
`cpdna_population_sites` track and mtDNA to the `mtdna_high_confidence_unique`
track. Calls are variant-only (`bcftools call -m -v`). Filtering and consensus
generation happen in later steps.

Run the controlled smoke call first:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_variant_calling.py \
  --run-label smoke \
  --sample-id ABAB_MAD_LP_223_Du-590 \
  --sample-id ABAB_MAD_LP_322_Du-593 \
  --sample-id ABAB_MAD_LP_323_Du-594 \
  --sample-id ABAB_MAD_LP_324_Du-595 \
  --sample-id ABAB_MAD_LP_325_Du-596
```

## Stage 19: Bootstrap-Supported Phylogenetic Trees

Stage 19 builds cpDNA and mtDNA maximum-likelihood trees from the full
callable-site consensus alignments. The final deliverable run uses IQ-TREE
with 1,000 ultrafast bootstrap replicates and BNNI correction.

Run from the repository root:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  MPLCONFIGDIR=/tmp/dudleya_matplotlib \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_phylogenetic_tree.py \
  --run-label primary \
  --output-dir dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree \
  --threads 4 \
  --bootstrap-replicates 1000
```

Default final tree outputs:

```text
dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile
dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile
dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/primary.phylogenetic_tree_summary.tsv
dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/primary.phylogenetic_tree_report.md
```

## Stage 15: PCA Visualization

Stage 15 computes cpDNA and mtDNA PCA from the filtered haploid SNP-only
alignments. It writes coordinate tables, variance summaries, and PNG/PDF/SVG
figures for each organelle.

Run from the repository root:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  MPLCONFIGDIR=/tmp/dudleya_matplotlib \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_pca_analysis.py \
  --run-label primary
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.coordinates.tsv
dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.variance.tsv
dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.png
dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.pdf
dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.svg
dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.coordinates.tsv
dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.variance.tsv
dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.png
dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.pdf
dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.svg
dudleya_organelle_alignment_pipeline/results/15_pca/primary.pca_summary.tsv
dudleya_organelle_alignment_pipeline/results/15_pca/primary.pca_report.md
```

## Stage 18: Admixture-Style Clustering

Stage 18 runs ADMIXTURE separately for cpDNA and mtDNA across a fixed K range
with cross-validation. For the final stability run, use multiple seeded
replicates per K and select the lowest mean CV error. Haploid organelle calls
are encoded as pseudo-diploid homozygotes for the diploid-oriented ADMIXTURE
tool, so the plots should be interpreted as organelle haplotype clustering
rather than nuclear admixture.

Run from the repository root:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  MPLCONFIGDIR=/tmp/dudleya_matplotlib \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_admixture_analysis.py \
  --run-label primary \
  --output-dir dudleya_organelle_alignment_pipeline/results/18_admixture_replicates \
  --max-k 8 \
  --threads 4 \
  --replicates 5
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/primary.admixture_summary.tsv
dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/primary.admixture_report.md
dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png
dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png
```

## Stage 14: Tree Visualization

Stage 14 renders the Stage 12 IQ-TREE Newick trees into static figures for
inspection and reporting. It does not alter the inferred tree topology. The final
publication figures are the Stage 20 renderings of the Stage 19
bootstrap-supported trees.

Run from the repository root:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  MPLCONFIGDIR=/tmp/dudleya_matplotlib \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_tree_visualization.py \
  --run-label primary
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/14_tree_visualization/cpDNA.primary.iqtree_ml_tree.png
dudleya_organelle_alignment_pipeline/results/14_tree_visualization/cpDNA.primary.iqtree_ml_tree.pdf
dudleya_organelle_alignment_pipeline/results/14_tree_visualization/cpDNA.primary.iqtree_ml_tree.svg
dudleya_organelle_alignment_pipeline/results/14_tree_visualization/mtDNA.primary.iqtree_ml_tree.png
dudleya_organelle_alignment_pipeline/results/14_tree_visualization/mtDNA.primary.iqtree_ml_tree.pdf
dudleya_organelle_alignment_pipeline/results/14_tree_visualization/mtDNA.primary.iqtree_ml_tree.svg
dudleya_organelle_alignment_pipeline/results/14_tree_visualization/primary.tree_visualization_summary.tsv
dudleya_organelle_alignment_pipeline/results/14_tree_visualization/primary.tree_visualization_report.md
```

## Additive R Visualization Alternatives

The pipeline also provides R-rendered alternatives for every existing PCA,
ADMIXTURE structure, ADMIXTURE cross-validation, initial-tree, and
bootstrap-tree figure. These files are additions: they read existing TSV and
Newick outputs, do not rerun the underlying biological analyses, and do not
replace or overwrite the original Matplotlib figures.

Run all 14 R figure jobs from the repository root:

```bash
python3 dudleya_organelle_alignment_pipeline/scripts/run_r_visualizations.py \
  --rscript .tools/bioconda-env/bin/Rscript
```

Pass `--force` to regenerate R alternatives that already exist, or select one
or more source stages with `--stages`, for example:

```bash
python3 dudleya_organelle_alignment_pipeline/scripts/run_r_visualizations.py \
  --rscript .tools/bioconda-env/bin/Rscript \
  --stages 15_pca 20_bootstrap_tree_visualization \
  --force
```

The additional filenames make the rendering method explicit:

```text
*.pca.r_ggplot.{png,pdf,svg}
*.structure.r_ggplot.{png,pdf,svg}
*.admixture_cv.r_ggplot.{png,pdf,svg}
*.iqtree_ml_tree.r_ggtree.{png,pdf,svg}
```

Legend and interpretation rules are generated with each graph:

- PCA and tree colors use the same fixed, colorblind-friendly species palette;
  unresolved metadata is gray. PCA no longer cycles 36 population groups
  through a 20-color palette.
- ADMIXTURE structure plots label colors as inferred clusters and state that
  cluster numbers and colors are arbitrary, not named biological populations.
  White boundaries and x-axis labels show the metadata population ordering.
- Cross-validation plots state that lower error is better, mark the selected K,
  and show mean plus/minus one standard deviation when replicates are present.
- Bootstrap tree figures identify internal numbers as UFBoot support
  percentages from 1,000 replicates. Tree tip colors are metadata annotations
  and do not affect inference.

Each affected results directory receives
`primary.r_visualization_commands.tsv` and
`primary.r_visualization_report.md` so the exact inputs, commands, outputs,
legend conventions, and interpretation limits remain auditable. The R
renderers use `ggplot2`, `ape`, and `ggtree`; the Stage 21 network renderer
continues to use `ape` and `pegas`.

## Stage 17: Fst And Population Summaries

Stage 17 computes pairwise population Fst and per-population summary statistics
from the filtered haploid SNP alignments. It uses only samples with resolved
population codes for population-level summaries.

Run from the repository root:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_population_genetics.py \
  --run-label primary
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv
dudleya_organelle_alignment_pipeline/results/17_population_genetics/cpDNA.primary.population_genetics.population_summary.tsv
dudleya_organelle_alignment_pipeline/results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv
dudleya_organelle_alignment_pipeline/results/17_population_genetics/mtDNA.primary.population_genetics.population_summary.tsv
dudleya_organelle_alignment_pipeline/results/17_population_genetics/primary.population_genetics_summary.tsv
dudleya_organelle_alignment_pipeline/results/17_population_genetics/primary.population_genetics_report.md
```

## Stage 21: Haploid Haplotype Networks

Stage 21 adds a haploid-native view of organelle sequence relationships without
removing or replacing the Stage 18 ADMIXTURE-style analysis. It reads the
filtered Stage 10 SNP alignments and Stage 07 sample metadata, then processes
cpDNA and mtDNA separately with `ape` and `pegas::haploNet`.

Before network inference, conservative complete-case filtering removes every
SNP column that contains a non-ACGT state in any sample. This keeps all 275
samples but drops sites with missing calls: cpDNA retains 1,977 of 2,015 sites
and mtDNA retains 116 of 146. The resulting networks contain 151 cpDNA
haplotypes and 58 mtDNA haplotypes.

Run from the repository root with the pipeline environment active:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  R_DEFAULT_DEVICE=png \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_haplotype_network.py \
  --run-label primary
```

In each figure, node area is proportional to the number of samples assigned to
the haplotype and colored sectors show species-group composition; unresolved
species metadata is gray. Exact sample-to-haplotype assignments, coordinates,
mutation counts, and primary versus alternative link status are retained in
TSV files. To prevent thousands of valid alternative connections from turning
the cpDNA panel into an unreadable hairball, the figures show primary links
only. For networks with more than 25 nodes, node IDs are omitted and primary
edges at or above the 90th percentile in mutation count (with a minimum of five
steps) are labeled. These are display rules only and do not change `haploNet`
inference.

Key outputs:

```text
dudleya_organelle_alignment_pipeline/results/21_haplotype_network/cpDNA.primary.haplotype_network.png
dudleya_organelle_alignment_pipeline/results/21_haplotype_network/mtDNA.primary.haplotype_network.png
dudleya_organelle_alignment_pipeline/results/21_haplotype_network/cpDNA.primary.haplotype_assignments.tsv
dudleya_organelle_alignment_pipeline/results/21_haplotype_network/mtDNA.primary.haplotype_assignments.tsv
dudleya_organelle_alignment_pipeline/results/21_haplotype_network/cpDNA.primary.haplotype_network_edges.tsv
dudleya_organelle_alignment_pipeline/results/21_haplotype_network/mtDNA.primary.haplotype_network_edges.tsv
dudleya_organelle_alignment_pipeline/results/21_haplotype_network/cpDNA.primary.popart.nex
dudleya_organelle_alignment_pipeline/results/21_haplotype_network/mtDNA.primary.popart.nex
dudleya_organelle_alignment_pipeline/results/21_haplotype_network/primary.haplotype_network_report.md
```

The PopART-compatible NEXUS files contain the filtered DNA matrix and one-hot
species trait columns for interoperability; the checked-in PNG/PDF/SVG figures
are rendered reproducibly by `pegas`. Because each organelle is inherited as a
linked lineage, these networks describe haplotype identity and sequence
connections. They are not ancestry proportions and do not by themselves infer
gene-flow direction, admixture fractions, divergence times, or nuclear
population structure.

Full-run command, after the smoke output is reviewed:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_variant_calling.py
```

Default outputs:

```text
dudleya_organelle_alignment_pipeline/results/08_variant_calling/cpDNA.raw.vcf.gz
dudleya_organelle_alignment_pipeline/results/08_variant_calling/mtDNA.raw.vcf.gz
dudleya_organelle_alignment_pipeline/results/08_variant_calling/variant_calling_summary.tsv
dudleya_organelle_alignment_pipeline/results/08_variant_calling/variant_calling_report.md
dudleya_organelle_alignment_pipeline/results/08_variant_calling/commands.tsv
```

Useful Stage 02 controlled-run options:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_pilot_alignment.py \
  --sample-limit 3

env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_pilot_alignment.py \
  --sample-id BAI_LP_105_Du-222

env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_pilot_alignment.py \
  --refresh-qc
```

To reproduce Stages 00 through 08:

```bash
python3 -m unittest \
  dudleya_organelle_alignment_pipeline.tests.test_manifest \
  dudleya_organelle_alignment_pipeline.tests.test_prepare_reference_and_pilot \
  dudleya_organelle_alignment_pipeline.tests.test_pilot_alignment \
  dudleya_organelle_alignment_pipeline.tests.test_analysis_masks \
  dudleya_organelle_alignment_pipeline.tests.test_all_sample_alignment \
  dudleya_organelle_alignment_pipeline.tests.test_downstream_sample_set \
  dudleya_organelle_alignment_pipeline.tests.test_variant_calling -v
python3 dudleya_organelle_alignment_pipeline/scripts/build_sample_manifest.py
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/prepare_reference_and_pilot.py
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_pilot_alignment.py
python3 dudleya_organelle_alignment_pipeline/scripts/build_analysis_masks.py
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_all_sample_alignment.py
python3 dudleya_organelle_alignment_pipeline/scripts/build_downstream_sample_set.py
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_variant_calling.py \
  --run-label smoke \
  --sample-id ABAB_MAD_LP_223_Du-590 \
  --sample-id ABAB_MAD_LP_322_Du-593 \
  --sample-id ABAB_MAD_LP_323_Du-594 \
  --sample-id ABAB_MAD_LP_324_Du-595 \
  --sample-id ABAB_MAD_LP_325_Du-596
```
