# Haploid Organelle Haplotype Networks Design

**Date:** 2026-07-14

## Objective

Add a reproducible Stage 21 that generates cpDNA and mtDNA haplotype-network
tables and figures with the established R population-genetics package `pegas`.
This is an additive analysis: no ADMIXTURE code, output, documentation, or final
deliverable is removed.

The new analysis addresses a different question from ADMIXTURE. It visualizes
how observed haploid organelle sequence types are connected by mutational
distance; it does not estimate nuclear ancestry, admixture proportions, or a
biological number of populations.

## Scientific Method

Stage 21 reads the filtered haploid SNP alignments produced by Stage 10 and
processes cpDNA and mtDNA separately. Before defining haplotypes, it removes any
alignment column containing a non-ACGT state in any sample. This complete-case
site filter keeps all 275 samples while preventing missing-data patterns from
being mistaken for biological haplotypes or silently imputed.

Based on the checked-in primary alignments, the expected inputs are:

| Organelle | Input SNP sites | Complete-case sites | Samples | Expected haplotypes |
|---|---:|---:|---:|---:|
| cpDNA | 2,015 | 1,977 | 275 | 151 |
| mtDNA | 146 | 116 | 275 | 58 |

The R rendering step uses `ape` to read DNA and `pegas` to collapse identical
sequences and construct a haplotype network. Node area represents the number of
samples carrying a haplotype. Node sectors show the six species/subspecies
metadata groups, including an `unresolved` group. Population-code composition is
preserved in machine-readable tables rather than encoded as 35 figure colors,
which would make the main diagrams unreadable.

Edges represent connections inferred by `pegas::haploNet` under its default
infinite-sites/Hamming-distance framework. Mutation steps are displayed using
the package's mutation-marking convention. Network layout is generated with a
fixed seed, and plotting coordinates are recorded in the output table so reruns
are reproducible. If `haploNet` cannot integrate the observed divergence, the
stage fails explicitly; it does not silently substitute a different graph method.

The report will state that:

- each organelle is a separate, largely linked cytoplasmic locus;
- network nodes are observed haplotypes after complete-case site filtering;
- edges are `haploNet` connections under the stated distance/model settings, not
  observed ancestor-descendant events;
- node sectors summarize sample metadata, not ancestry proportions;
- cpDNA and mtDNA networks must be interpreted separately;
- the networks supplement, rather than replace, PCA, trees, Fst, and the retained
  exploratory ADMIXTURE outputs.

## Architecture

### Python stage module

Create `dudleya_organelle_alignment_pipeline/haplotype_network.py`. It follows
the existing stage-module pattern and owns:

- reading the Stage 10 summary and SNP FASTA files;
- reading the Stage 07 included-sample metadata;
- validating sample identity and order;
- complete-case site filtering;
- writing filtered network-input FASTA files;
- writing sample metadata for R;
- building and running the `Rscript` command;
- validating required R output files;
- writing the combined Stage 21 summary and Markdown report.

Pure data preparation remains in Python so it is covered by the repository's
existing `unittest` suite and shares the same input-validation conventions as
other stages.

### R network renderer

Create `dudleya_organelle_alignment_pipeline/scripts/render_haplotype_network.R`.
It owns only domain-specific network calculation and rendering:

- read the complete-case FASTA with `ape`;
- build haplotypes and a network with `pegas`;
- join sample metadata to haplotypes;
- write assignment, haplotype-frequency, edge, and layout tables;
- render PNG, PDF, and SVG figures with species-composition node sectors.

The script exits nonzero with a specific message when required packages, input
columns, samples, sites, or graph outputs are missing.

### Runner

Create `dudleya_organelle_alignment_pipeline/scripts/run_haplotype_network.py`.
It exposes the normal pipeline arguments: input summary, sample metadata,
output directory, run label, Rscript executable, renderer path, and force flag.
Defaults point to the checked-in primary Stage 10 and Stage 07 outputs and write
to `results/21_haplotype_network/`.

## Stage 21 Outputs

For each organelle, Stage 21 writes:

- `<organelle>.primary.haplotype_network_input.fa`
- `<organelle>.primary.haplotype_network_sites.tsv`
- `<organelle>.primary.haplotype_assignments.tsv`
- `<organelle>.primary.haplotype_summary.tsv`
- `<organelle>.primary.haplotype_network_edges.tsv`
- `<organelle>.primary.haplotype_network_layout.tsv`
- `<organelle>.primary.haplotype_network.png`
- `<organelle>.primary.haplotype_network.pdf`
- `<organelle>.primary.haplotype_network.svg`
- `<organelle>.primary.popart.nex`

The stage also writes:

- `primary.haplotype_network_summary.tsv`
- `primary.haplotype_network_commands.tsv`
- `primary.haplotype_network_report.md`

The PopART-compatible NEXUS files are supplementary interoperability products.
The checked-in figures are generated by `pegas`, not manually through PopART.

## Dependencies And Tool Audit

Add `r-pegas` to `environment.yml` and add an `r_pegas` import check to the Stage
13 tool audit. `Rscript`, `ape`, and the existing Python plotting stack remain
required. Documentation records the observed `pegas` version used for the final
run.

## Documentation

Update only pipeline and result documentation, not the study guide:

- `dudleya_organelle_alignment_pipeline/PROCESS.md` adds Stage 21 and changes the
  documented contiguous stage range to 00-21.
- `dudleya_organelle_alignment_pipeline/README.md` documents purpose, method,
  command, outputs, interpretation, and reproduction.
- `dudleya_organelle_alignment_pipeline/results/organelle_population_report.md`
  adds the haplotype-network methods, results, figures, and limitations while
  retaining the ADMIXTURE section.
- `dudleya_organelle_alignment_pipeline/results/final_deliverables_manifest.tsv`
  adds both network figures, both assignment tables, and the Stage 21 report.
- `dudleya_organelle_alignment_pipeline/results/13_tool_audit/` is regenerated
  after the dependency is installed.

## Error Handling

The stage fails before R execution if input files are missing, FASTA lengths are
inconsistent, samples differ from the included-sample table, no complete-case
sites remain, or fewer than two haplotypes remain. It fails after R execution if
the command exits nonzero or any required table or figure is absent or empty.

Temporary figure files are written before final replacement so an interrupted
run cannot leave an apparently complete deliverable. Exact R commands are stored
in the Stage 21 command table.

## Testing And Verification

Development follows test-first cycles. Add
`tests/test_haplotype_network.py` covering:

- complete-case filtering removes every column with `N` or another non-ACGT
  state while retaining all samples;
- retained-site coordinates are recorded correctly;
- sample order/identity mismatches fail loudly;
- R command construction includes all required paths and run labels;
- required-output validation detects missing and empty products;
- summary and report writers record sample, site, haplotype, and figure counts.

Add a tool-audit test proving `r_pegas` is checked. After unit tests pass, run the
primary Stage 21 analysis, inspect both diagrams, verify all tables against the
source alignments, rerun the full test suite, and check documentation links and
the final-deliverables manifest.

## Non-Goals

- Removing, renaming, or recalculating existing ADMIXTURE outputs.
- Interpreting network sectors as ancestry or mixture proportions.
- Combining cpDNA and mtDNA into one network.
- Imputing missing SNP calls.
- Inferring ancestral haplotypes, migration direction, species boundaries, or a
  species tree from the network alone.
- Editing any file in `dudleya_organelle_alignment_pipeline/study_guide/`.
