# Chapter 18 — Exercises

> Part 4 of 4 · Practice and Reference · Prev:
> [Uncertainty, Bias, and Limits](./17-uncertainty-bias-and-limits.md) · Next:
> [Solutions](./19-solutions.md)

Two kinds of exercises. Each Part 2 walkthrough chapter ends with six exercises
tied to that stage; this chapter adds a warm-up, a set of *integrative* exercises
that cross stage boundaries, and a guide to learning from the test suite. All
solutions are in [Chapter 19](./19-solutions.md), keyed by number.

The five exercise types recur throughout, matching the chapter teaching pattern:

- **Trace** — run code in your head on given input.
- **Predict** — say what an output will be before checking.
- **Modify** — make a small, correct code or test change.
- **Debug** — diagnose a described failure.
- **Interpret** — audit a biological or statistical claim.

## 18.0 Warm-up: run the tests

Before anything else, from the repository root:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 -m unittest discover -s dudleya_organelle_alignment_pipeline/tests -v
```

**W1.** How many test files and how many tests run, and roughly how long does the
suite take? Why does it not need `bwa`, `iqtree`, or `admixture` installed — and
which tests still need non-standard Python libraries in the active environment?

**W2.** Run just one module's tests:
`python3 -m unittest dudleya_organelle_alignment_pipeline.tests.test_callable_consensus -v`.
Which test asserts the exact strings `"CTTNNG"` and `"CGTANG"`, and what stage
does it pin?

## 18.1 Per-chapter exercises (index)

Each stage chapter has its own six-exercise set:

- Stages 00–01: [Chapter 7, §7.9](./07-manifest-and-reference-preflight.md)
- Stages 02–04: [Chapter 8, §8.8](./08-pilot-mapping-and-investigations.md)
- Stages 05–07: [Chapter 9, §9.7](./09-masks-alignment-and-sample-qc.md)
- Stages 08–11: [Chapter 10, §10.9](./10-variants-to-alignments.md)
- Stages 12/14/19/20: [Chapter 11, §11.8](./11-phylogenetic-trees.md)
- Stages 15/16/18: [Chapter 12, §12.6](./12-pca-and-clustering.md)
- Stage 17: [Chapter 13, §13.9](./13-population-fst.md)
- Stage 13: [Chapter 14, §14.9](./14-tool-audit.md)

## 18.2 Integrative exercises (cross-stage)

These require connecting two or more stages.

**I1. Trace a number through the chain.** Starting from 2,475 raw cpDNA variant
records, name every stage that transforms that number and give the value after
each, ending at the SNP-alignment column count. Which single `bcftools view`
flag removes the most records, and why?

**I2. Predict a break.** You rerun Stage 09 with `--max-missing-fraction 0.5` and
`--min-minor-allele-count 1`. Directionally, what happens to the cpDNA filtered
SNP count, the Stage 10 alignment width, and the Stage 15 PCA `retained_sites`?
Would any stage *raise*?

**I3. The run label.** You run Stage 08 with `--run-label smoke` and five
samples, then run Stage 09 with the default `--run-label primary`. What goes
wrong, and which function's naming logic explains it?

**I4. Track discipline.** A colleague edits `analysis_tracks.tsv` so that variant
calling points at `cpdna_full_coverage` instead of `cpdna_population_sites`.
Nothing crashes. What biological artifact returns, roughly how would the cpDNA
SNP count change, and which chapter's rule was violated?

**I5. Modify across two stages.** You want a ≥20× breadth QC metric to appear in
both the pilot (Stage 02) and all-sample (Stage 06) reports. List every function
or data structure you must touch in each module. Which alignment machinery is
shared, and which track-depth machinery is Stage-06-specific?

**I6. Debug a coverage paradox.** A sample has `cpdna_full_coverage` breadth 0.95
but `cpdna_population_sites` breadth 0.99, and its callable-consensus record is
mostly `N`. Reconcile the two breadth numbers, then explain why “mostly `N`” is
inconsistent with 0.99 breadth at the same ≥1× threshold and what configuration
or file mismatch you would investigate.

**I7. Interpret a full result.** For cpDNA, populations P and Q sit far apart on
PC1, fall in separate UFBoot-98 clades, and show pairwise Fst 0.45; for mtDNA
they overlap on PC1, sit in a low-support region of the mtDNA tree, and show Fst
0.02. Write the single most defensible sentence describing what this supports,
using the evidence tags and the one-locus caveat.

**I8. Capstone rehearsal.** Using the six-base example from the callable-consensus
test (reference `ACGTACGT`, BED `chloroplast 1 7`, one filtered SNP at position
3, one raw-only failed site at position 6, and the two samples' depth files),
hand-compute sample S1's consensus string and confirm it against the test.
(The full trace is [Chapter 23](./23-capstone-sample-trace.md) — try it first.)

## 18.3 Learning from the tests directly

The fastest way to deepen your understanding is to read a test, predict what it
asserts, then run it. Try these:

**T1.** In `test_analysis_masks.py`, the cpDNA fixture uses IR copies at `5–8`
and `15–18` on a 20-bp reference. Before reading the assertion, predict the
`cpdna_population_sites` regions. Then confirm.

**T2.** In `test_variant_filtering.py`, list every `bcftools view` flag the
command-builder test asserts is present, and match each to a biological or
quality reason from [Chapter 10](./10-variants-to-alignments.md).

**T3.** In `test_admixture_analysis.py`, `summarize_replicate_stability` is given
K=2 (mean 0.31) and K=3 (mean 0.21). Predict which K is `is_best_mean_k`, then
change the fixture so K=2 wins and predict the new assertion.

**T4.** Pick any module and add one new test: a small fixture, one function call,
one assertion. State what behavior you are pinning and why it was not already
covered.

## 18.4 A modification project (optional, larger)

If you want a bigger exercise: add a `--min-depth`-configurable *sample-level*
callability summary that reports, per sample, the fraction of its
callable-consensus that is non-`N`. Sketch which module it belongs in
([`callable_consensus.py`](../callable_consensus.py)), which existing function
already computes the missing count you would reuse, the new TSV columns, and the
test you would write first (test-driven: write the failing assertion, then the
code). You do not need to run it — the design and the test are the exercise.

> Next: [Chapter 19 — Solutions](./19-solutions.md)
