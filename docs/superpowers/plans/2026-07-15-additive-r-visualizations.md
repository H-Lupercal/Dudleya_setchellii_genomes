# Additive R Visualizations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate additive R versions of every existing non-R analysis figure with automatic legends or explanatory keys.

**Architecture:** A Python orchestration module discovers existing analysis outputs and invokes one R renderer per figure family. The R scripts read existing TSV/Newick data, use a shared species palette, and write renderer-suffixed PNG, PDF, and SVG files without touching current figures.

**Tech Stack:** Python 3, R 4.5, ggplot2, ggtree, ape, patchwork, unittest.

---

### Task 1: Define orchestration behavior with failing tests

**Files:**
- Create: `dudleya_organelle_alignment_pipeline/tests/test_r_visualizations.py`
- Create: `dudleya_organelle_alignment_pipeline/r_visualizations.py`

- [ ] **Step 1: Write tests for additive paths and renderer commands**

Create tests that import `figure_outputs`, `build_renderer_command`, and
`discover_figure_jobs`. Assert that PCA, structure, CV, and tree jobs end in
`.r_ggplot` or `.r_ggtree`, contain PNG/PDF/SVG paths, and pass existing source
files rather than existing image paths to R.

- [ ] **Step 2: Run the new tests and verify RED**

Run:

```bash
python3 -m unittest dudleya_organelle_alignment_pipeline.tests.test_r_visualizations -v
```

Expected: import failure because `r_visualizations.py` does not exist.

- [ ] **Step 3: Implement the minimal orchestration API**

Add dataclasses for figure jobs and results; implement deterministic output
naming, input discovery for Stages 14, 15, 16, 18, and 20, R command creation,
subprocess execution, and command/report TSV/Markdown writers.

- [ ] **Step 4: Run the new tests and verify GREEN**

Run the Task 1 unittest command. Expected: all new orchestration tests pass.

### Task 2: Add shared plotting rules and PCA renderer test-first

**Files:**
- Create: `dudleya_organelle_alignment_pipeline/scripts/dudleya_plotting.R`
- Create: `dudleya_organelle_alignment_pipeline/scripts/render_pca_ggplot.R`
- Modify: `dudleya_organelle_alignment_pipeline/tests/test_r_visualizations.py`

- [ ] **Step 1: Add a failing R smoke test**

Write a temporary coordinate TSV with `sample_id`, `pc1`, `pc2`, `species`, and
`plot_group`, plus a two-row variance TSV. Run the intended renderer and assert
that all three additive output files are non-empty.

- [ ] **Step 2: Verify the smoke test fails because the renderer is missing**

Run the targeted unittest and confirm the missing-script failure.

- [ ] **Step 3: Implement shared palette/output helpers and PCA plotting**

Use the six species categories, convert blanks to `unresolved`, label PC axes
from the variance TSV, and use a right-side `Species group` legend.

- [ ] **Step 4: Verify the PCA smoke test passes**

Run the targeted unittest and confirm PNG/PDF/SVG outputs are non-empty.

### Task 3: Add ADMIXTURE structure and CV renderers test-first

**Files:**
- Create: `dudleya_organelle_alignment_pipeline/scripts/render_admixture_ggplot.R`
- Modify: `dudleya_organelle_alignment_pipeline/tests/test_r_visualizations.py`

- [ ] **Step 1: Add failing structure and CV smoke tests**

Provide a small ordered Q TSV and summary TSV. Assert that structure and CV
modes each write three files and reject missing cluster columns.

- [ ] **Step 2: Verify RED**

Run the targeted tests and confirm failure because the renderer is absent.

- [ ] **Step 3: Implement both modes**

Structure mode pivots `cluster_*` columns into stacked assignment bars, draws
population boundaries, and titles the fill legend `Inferred cluster\n(labels arbitrary)`.
CV mode collapses replicates by K, draws mean plus standard deviation, marks the
minimum-mean K, and adds `Lower cross-validation error is better`.

- [ ] **Step 4: Verify GREEN**

Run the targeted tests and confirm both modes pass.

### Task 4: Add ggtree renderer test-first

**Files:**
- Create: `dudleya_organelle_alignment_pipeline/scripts/render_tree_ggtree.R`
- Modify: `dudleya_organelle_alignment_pipeline/tests/test_r_visualizations.py`

- [ ] **Step 1: Add failing initial-tree and bootstrap-tree smoke tests**

Create a small Newick tree and matching metadata. Assert both modes generate
three formats and that bootstrap mode accepts a replicate count.

- [ ] **Step 2: Verify RED**

Run the targeted tests and confirm the missing renderer failure.

- [ ] **Step 3: Implement ggtree plotting**

Join tip labels to metadata, color tip labels by species group, preserve branch
lengths, show support labels in bootstrap mode, and add the UFBoot explanation.

- [ ] **Step 4: Verify GREEN**

Run the targeted tests and confirm both modes pass.

### Task 5: Add CLI, generate artifacts, and document them

**Files:**
- Create: `dudleya_organelle_alignment_pipeline/scripts/run_r_visualizations.py`
- Modify: `dudleya_organelle_alignment_pipeline/README.md`
- Create: additive R figures and `primary.r_visualization_{commands.tsv,report.md}` in the existing results directories

- [ ] **Step 1: Add a failing CLI/parser test**

Assert defaults point at the repository result and script directories and that
`--rscript`, `--run-label`, and selected stage flags are accepted.

- [ ] **Step 2: Verify RED, then implement the CLI and verify GREEN**

The CLI calls `run_r_visualizations`, prints each added artifact, and never
deletes or overwrites an existing non-R figure.

- [ ] **Step 3: Generate all additive figures**

Run:

```bash
python3 dudleya_organelle_alignment_pipeline/scripts/run_r_visualizations.py --rscript .tools/bioconda-env/bin/Rscript
```

Expected: R alternatives for 2 PCA, 8 ADMIXTURE/CV, and 4 tree figures, each in
PNG/PDF/SVG.

- [ ] **Step 4: Document commands, filenames, legends, and limits**

Add a README section explaining that these are visualization alternatives, not
new biological analyses, and that existing figures remain canonical inputs.

### Task 6: Verify the complete additive change

**Files:**
- Verify all files above and every new PNG.

- [ ] **Step 1: Run targeted and full tests**

```bash
python3 -m unittest dudleya_organelle_alignment_pipeline.tests.test_r_visualizations -v
python3 -m unittest discover -s dudleya_organelle_alignment_pipeline/tests -v
```

Expected: zero failures.

- [ ] **Step 2: Verify artifact inventory and original preservation**

Confirm all expected additive files are non-empty, no original figure is
modified, and every new basename contains `.r_ggplot` or `.r_ggtree`.

- [ ] **Step 3: Visually inspect every new PNG**

Check titles, legends/keys, axis labels, support explanations, population
boundaries, clipping, and readability.

- [ ] **Step 4: Run repository hygiene checks**

```bash
git diff --check
git status --short
```

Expected: no whitespace errors and additions/modifications limited to this
feature. Do not commit or push.
