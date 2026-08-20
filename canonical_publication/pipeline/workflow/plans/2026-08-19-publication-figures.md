# Canonical Publication Figures Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a fingerprinted, deterministic publication-figure stage that renders the approved organelle figure suite only from canonical outputs.

**Architecture:** Analysis stages retain responsibility for data products; a new renderer module owns tested visual primitives and a new stage script composes publication panels. The workflow inserts this stage between supplementary ADMIXTURE and reports, and report acceptance validates its state and manifest.

**Tech Stack:** Python 3.14, Matplotlib, NumPy, pandas, NetworkX, Biopython, pytest, Ruff, canonical provenance helpers.

---

### Task 1: Lock visual semantics in pure helpers

**Files:**
- Create: `canonical_publication/pipeline/src/organelle_pipeline/figures.py`
- Create: `canonical_publication/pipeline/tests/test_figures.py`

- [ ] **Step 1: Write failing tests** for the five-taxon palette, symmetric signed-FST color limits around zero, deterministic distance-aware haplotype layout, taxon pie counts, and unrooted tree coordinates.
- [ ] **Step 2: Verify red** with `pytest -q canonical_publication/pipeline/tests/test_figures.py`; expect import failure because `organelle_pipeline.figures` does not exist.
- [ ] **Step 3: Implement the smallest pure helpers** that validate inputs and return deterministic colors, bounds, counts, and coordinates without writing files.
- [ ] **Step 4: Verify green** with `pytest -q canonical_publication/pipeline/tests/test_figures.py`; expect all figure-helper tests to pass.
- [ ] **Step 5: Re-run the full unit suite** with `pytest -q canonical_publication/pipeline/tests`; expect no regressions.

### Task 2: Separate analysis data from presentation

**Files:**
- Modify: `canonical_publication/pipeline/scripts/run_pca.py`
- Modify: `canonical_publication/pipeline/scripts/run_haplotypes.py`
- Modify: `canonical_publication/pipeline/tests/test_analysis_commands.py`

- [ ] **Step 1: Add failing source/behavior tests** asserting that PCA and haplotype states declare only tabular outputs and no longer create Matplotlib figures.
- [ ] **Step 2: Verify red** with the targeted pytest nodes; expect assertions to find current PNG/PDF/SVG handling.
- [ ] **Step 3: Remove renderer imports, directories, figure paths, plotting, and figure checksums** from the two analysis scripts while preserving numerical outputs and fingerprints.
- [ ] **Step 4: Verify green** with the targeted tests and then the full suite.

### Task 3: Implement the fingerprinted renderer

**Files:**
- Create: `canonical_publication/pipeline/scripts/render_figures.py`
- Extend: `canonical_publication/pipeline/tests/test_figures.py`

- [ ] **Step 1: Add failing miniature integration tests** that create canonical fixture tables/trees/states and require nonempty `.png`, `.pdf`, `.svg`, a `figure_manifest.tsv`, and a provenance JSON whose outputs match file SHA-256 values.
- [ ] **Step 2: Verify red**; expect failure because `render_figures.py` does not exist.
- [ ] **Step 3: Implement input-state validation and panel renderers** for reference/callability, preprocessing/QC, PCA, haplotype networks, signed Hudson FST, unrooted primary trees, and supplementary ADMIXTURE.
- [ ] **Step 4: Implement atomic figure-set publication** to the run-specific figure directory, write the manifest, and save a stage fingerprint containing all upstream state digests and renderer/runtime versions.
- [ ] **Step 5: Verify green** with the miniature integration test and full test suite.

### Task 4: Integrate figures into dependency ordering and acceptance

**Files:**
- Modify: `canonical_publication/pipeline/scripts/run_pipeline.py`
- Modify: `canonical_publication/pipeline/scripts/build_reports.py`
- Modify: `canonical_publication/pipeline/tests/test_layout.py`
- Modify: `canonical_publication/pipeline/tests/test_workflow_commands.py`

- [ ] **Step 1: Add failing workflow tests** requiring `admixture < figures < reports`, the exact renderer command, report fingerprint inclusion of `figures.json`, and acceptance rejection when the figure manifest/state is absent.
- [ ] **Step 2: Verify red** with the targeted workflow/layout tests.
- [ ] **Step 3: Add the `figures` command and stage order**, then require figure outputs and provenance in report construction and acceptance.
- [ ] **Step 4: Verify green** with targeted tests and the full unit suite.

### Task 5: Static verification before expensive stages

**Files:**
- Modify only files implicated by failures from the commands below.

- [ ] **Step 1: Run Ruff checks** with `ruff check canonical_publication/pipeline` and `ruff format --check canonical_publication/pipeline`.
- [ ] **Step 2: Compile Python** with `python -m compileall -q canonical_publication/pipeline`.
- [ ] **Step 3: Parse shell** with `bash -n canonical_publication/run_pipeline.sh`.
- [ ] **Step 4: Check whitespace** with `git diff --check`.
- [ ] **Step 5: Dry-run the complete workflow** and confirm the displayed order includes `STAGE figures` immediately before `STAGE reports`.

### Task 6: Regenerate the canonical run and inspect figures

**Files:**
- Generate: `canonical_publication/reports/figures/publication-20260817/`
- Generate: `canonical_publication/provenance/runs/publication-20260817/figures.json`

- [ ] **Step 1: Run the complete resumable pipeline** with `SCIKIT_ALLEL_PYTHON=<validated-scikit-allel-python> canonical_publication/run_pipeline.sh --config canonical_publication/config/publication_config.toml --run-id publication-20260817 --resume`, supplying the validated interpreter through the environment rather than embedding a workstation path in a canonical artifact.
- [ ] **Step 2: Verify all upstream resume checks and downstream stages complete**, with no stale-output bypass.
- [ ] **Step 3: Inspect every generated PNG**, compare it to its PDF/SVG manifest entries, and correct only through a new test-first code change followed by fingerprint-aware regeneration.
- [ ] **Step 4: Confirm `ACCEPTANCE.json` reports PASS**, sample counts match QC-derived organelle tables, pair counts equal `n(n-1)/2`, and every final figure is checksummed.
