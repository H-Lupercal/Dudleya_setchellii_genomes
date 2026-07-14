# Dudleya Organelle Pipeline Study Guide Design

## Goal

Create a self-contained textbook and study guide for
`dudleya_organelle_alignment_pipeline` that teaches its Python, command-line
bioinformatics, and biology. The primary reader is strong in biology and has
some technical experience, but needs more support on Python and software
engineering. The guide must support both systematic reading and hands-on study.

## Scope

The guide covers the complete pipeline surface:

- Stages 00 through 20 and the initial-versus-final analysis distinction.
- Every Python module, runner script, and top-level function.
- The tests as executable statements of intended behavior.
- The checked-in reports, summaries, command logs, tables, alignments, trees,
  and figures needed to understand what the pipeline produced.
- The external tools and file formats used by the pipeline, at the depth needed
  to read and reason about the repository.
- The biological assumptions, inference limits, quality controls, and failure
  modes behind each stage.

The guide does not duplicate large FASTQ files, require a full pipeline rerun,
or claim that organelle analyses establish nuclear admixture or complete
species history.

## Reader Model

The reader has a strong biology background and wants pipeline-specific mastery.
Python concepts therefore receive explicit teaching and annotated examples,
while familiar biological vocabulary can be introduced more compactly. The
book still explains organelle-genomic concepts whose computational consequences
are easy to miss, including haploidy, uniparental inheritance, repeats,
callability, reference bias, linkage, and discordance between organelle trees.

## Deliverable Architecture

The textbook lives in
`dudleya_organelle_alignment_pipeline/study_guide/`. Its `README.md` is the book
home and table of contents. The pipeline's existing `README.md` receives a short
link to the book; its operational instructions remain authoritative and are not
rewritten.

The book has four parts:

1. **Foundations**: the full data-flow map, essential Python and shell concepts,
   bioinformatics file formats, and organelle biology.
2. **Pipeline walkthrough**: Stages 00-20, grouped into coherent biological and
   computational units rather than assigning every stage equal space.
3. **Interpretation**: phylogenetics, PCA, ADMIXTURE-style clustering, pairwise
   Fst, uncertainty, bias, and appropriate limits on biological claims.
4. **Practice and reference**: exercises, solutions, glossary, module/function
   index, external-tool reference, and a sample-tracing capstone.

All chapters are plain Markdown and use relative links to code, tests, runner
scripts, and representative outputs.

## Chapter Teaching Pattern

Each pipeline chapter follows the same sequence:

1. State the biological question or quality-control problem.
2. Identify the actual input and output files.
3. Map runner script to Python module, external commands, and results.
4. Walk through important functions and code blocks.
5. Teach the Python and software-engineering concepts used by that code.
6. Interpret the biological result and distinguish it from unsupported claims.
7. Explain biological and computational failure modes.
8. Provide code-reading, prediction, debugging, modification, and
   interpretation exercises.

Repeated helper patterns are taught once and cross-referenced. The reference
index still lists every module and top-level function so coverage remains
auditable.

## Data Flow

The guide presents the pipeline as a provenance chain:

```text
FASTQ files + population metadata
        |
        v
manifest and reference preflight (00-01)
        |
        v
pilot mapping and reference investigations (02-04)
        |
        v
masks, all-sample alignment, sample QC (05-07)
        |
        v
haploid calls, filters, SNP/callable alignments (08-11)
        |
        +----------------------+----------------------+------------------+
        v                      v                      v                  v
phylogenetic trees       PCA and clustering      population Fst     tool audit
(12, 14, 19, 20)         (15, 16, 18)            (17)               (13)
```

Initial tree and clustering runs are labeled as historical first passes. Final
interpretation points readers to the five-replicate clustering in Stage 18 and
the 1,000-replicate bootstrap trees in Stages 19-20.

## Code Coverage Strategy

The research pass reads every implementation module and test end-to-end. For
each module, the reference records:

- Its purpose and stage.
- Its public functions, parameters, return values, and raised errors.
- The files it reads and writes.
- External programs it invokes.
- Tests that define its expected and edge-case behavior.
- Biological assumptions embedded in thresholds, masks, encodings, or sample
  inclusion decisions.

The narrative focuses on consequential code paths. Small formatting and path
helpers are covered through a reusable-pattern chapter and the complete
function index.

## Exercises and Solutions

Exercises use small fixtures, existing unit tests, and checked-in summaries so
the reader can learn without rerunning the full 275-sample workflow. Each
chapter includes a mixture of:

- Trace-the-code questions.
- Predict-the-output questions.
- Small Python changes or test additions.
- Debugging scenarios based on validated failure paths.
- Biological interpretation and claim-auditing questions.

Solutions explain the reasoning and point back to source lines or tests. The
capstone follows one sample conceptually from paired FASTQs through QC, variant
representation, downstream matrices, and cautious biological interpretation.

## Error Handling and Scientific Caution

The guide explains both software failures and scientific failure modes:

- Missing tools, malformed paths, failed subprocesses, and partial outputs.
- Missing mates and invalid sample metadata.
- Malformed TSV, FASTA, BED, BAM, and VCF inputs.
- Low depth, low callability, excessive missingness, and sample exclusion.
- Repeats, paralogous placement, inverted repeats, and reference bias.
- Haploid calling assumptions and linked organelle markers.
- Overinterpretation of PCA, clustering, Fst, topology, and bootstrap support.

Claims in the book are labeled by evidence type: source-code behavior,
test-established behavior, observed repository result, or biological
interpretation.

## Verification

Before completion:

1. Run the existing unit-test suite.
2. Check every new relative Markdown link.
3. Compare the stage sequence against `PROCESS.md`.
4. Verify that every implementation module, runner, and top-level function is
   present in the reference index.
5. Check examples against actual signatures, defaults, and test fixtures.
6. Scan for placeholders, contradictions, unsupported biological claims, and
   ambiguous initial-versus-final result descriptions.

## Success Criteria

The finished guide lets the reader:

- Explain what each stage does, why it exists, and how data moves onward.
- Read every Python module without being blocked by unfamiliar syntax or
  architectural patterns.
- Understand how shell commands and external tools are assembled and audited.
- Connect filtering and QC choices to downstream biological inference.
- Interpret the repository's organelle results with appropriate caution.
- Use tests and small exercises to verify their understanding.
- Locate any module, function, tool, format, stage, or biological term quickly.

