# Additive R Visualizations Design

## Goal

Add publication-oriented R versions of the existing PCA, phylogenetic-tree,
ADMIXTURE structure, and ADMIXTURE cross-validation figures. Preserve every
existing script and result file.

## Approaches Considered

1. Replace the Python plots with R plots. This would create one canonical
   figure per analysis, but it violates the requirement to remove nothing and
   makes comparisons difficult.
2. Add one large R script that detects every input type. This minimizes the
   number of scripts, but couples unrelated plotting rules and makes failures
   harder to diagnose.
3. Add focused R renderers plus a small Python orchestrator. This keeps the
   biological analyses unchanged, gives each figure family an explicit input
   contract, and writes clearly named alternatives beside the originals.

Approach 3 is selected.

## Architecture and Data Flow

`r_visualizations.py` discovers existing result tables and tree files and calls
three focused renderers:

- `render_pca_ggplot.R` reads Stage 15 coordinate and variance TSV files.
- `render_admixture_ggplot.R` reads Stage 16 or 18 Q tables and ADMIXTURE
  summaries.
- `render_tree_ggtree.R` reads Stage 14 or 20 Newick trees and the Stage 7
  sample metadata.

The renderers do not recalculate PCA, clustering, likelihood trees, or
bootstrap values. They only visualize existing outputs. A command TSV and a
Markdown report document inputs, outputs, packages, legends, and biological
interpretation limits.

## Additive Output Names

Each new figure uses the original prefix plus a renderer suffix:

- PCA: `*.pca.r_ggplot.{png,pdf,svg}`
- structure: `*.structure.r_ggplot.{png,pdf,svg}`
- CV error: `*.admixture_cv.r_ggplot.{png,pdf,svg}`
- trees: `*.iqtree_ml_tree.r_ggtree.{png,pdf,svg}`

Existing filenames are never opened for writing.

## Visual and Scientific Rules

- Use one fixed, colorblind-friendly six-category species palette across PCA
  and trees, with blank species metadata displayed as `unresolved`.
- PCA color represents species group, not the 36 population codes. The legend
  is always present, and PC axes retain their explained-variance percentages.
- ADMIXTURE fills represent inferred clusters. The legend title states that
  cluster labels are arbitrary. Population boundaries are drawn and labeled
  from the existing ordered Q table.
- CV figures state that lower error is better, mark the selected K, and show
  mean plus standard deviation when replicate summaries are available.
- Bootstrap trees explain that internal values are UFBoot support percentages
  from 1,000 replicates. Initial trees omit that explanation because they do
  not contain bootstrap support values.
- Tree branch lengths remain in substitutions per site. Species colors annotate
  tips but do not change tree inference.
- Every new figure contains either a categorical legend or an explanatory
  on-figure key appropriate to its visual encodings.

## Error Handling

The orchestrator validates required inputs and R scripts before launching a
renderer. Each renderer fails on missing required columns, unknown figure mode,
or empty data. Subprocess failures include the command output and leave existing
figures untouched.

## Testing and Verification

- Unit tests cover input discovery, additive output naming, command creation,
  and report/manifest content.
- R smoke tests generate small PCA, ADMIXTURE, CV, and tree figures and verify
  all three formats are non-empty.
- The repository test suite and `git diff --check` run after generation.
- Every generated PNG is visually inspected for its title, legend/key,
  readable labels, and absence of clipping.

## Scope

The work adds visualization code, tests, documentation, and generated
artifacts only. It does not alter biological inputs, regenerate analyses,
replace existing figures, change study-guide chapters, commit, or push.
