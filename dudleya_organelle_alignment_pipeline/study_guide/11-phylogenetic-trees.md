# Chapter 11 — Phylogenetic Trees (Stages 12, 14, 19, 20)

> Part 2 of 4 · Pipeline Walkthrough · Prev:
> [Variants to Alignments](./10-variants-to-alignments.md) · Next: [PCA and
> Clustering](./12-pca-and-clustering.md)

Four stages, one module pair. [`phylogenetic_tree.py`](../phylogenetic_tree.py)
infers trees; [`tree_visualization.py`](../tree_visualization.py) renders them.
Stages 12/14 are the fast first pass; Stages 19/20 are the bootstrap-supported
final deliverable. This chapter is about *building* the trees; *reading* them is
[Chapter 15](./15-phylogenetics-interpretation.md).

## 11.1 The question

*Given the callable-site alignment for each organelle, what is the
maximum-likelihood tree of the samples, and how confident are we in each branch?*

## 11.2 The stages

| Stage | What | Runner arguments |
|---|---|---|
| 12 | ML trees, fast, no bootstrap | (defaults) |
| 14 | figures for Stage 12 trees | (defaults) |
| 19 | ML trees, 1,000 UFBoot + BNNI | `--bootstrap-replicates 1000 --output-dir .../19_...` |
| 20 | figures for Stage 19 trees | `--tree-dir .../19_...` |

Input: the callable-site consensus FASTAs from Stage 11 (275 × 124,538 cpDNA;
275 × 44,930 mtDNA). Output: IQ-TREE `.treefile` (Newick), a summary TSV, a
commands TSV, a `report.md`, and rendered PNG/PDF/SVG figures.

## 11.3 The code: one command, two modes

`build_iqtree_command` assembles the IQ-TREE invocation ([Chapter 4,
§4.6](./04-shell-and-external-tools.md)):

```python
command = [iqtree_executable, "-s", alignment_path.as_posix(),
           "--seqtype", "DNA", "-m", model, "--prefix", prefix.as_posix(),
           "-T", str(threads), "--safe", "--redo", "--quiet"]
if fast:
    command.append("--fast")
if bootstrap_replicates:
    command.extend(["-B", str(bootstrap_replicates), "--bnni"])
```

The model default is **`GTR+F+G4`**. The two modes are mutually exclusive, and
the switch lives in `main`:

```python
fast = not args.full_search and not args.bootstrap_replicates
```

So asking for bootstraps (`--bootstrap-replicates 1000`) automatically sets
`fast=False`, which is exactly what separates Stage 19 from Stage 12. The method
label recorded in the summary reflects this: `iqtree_ml_fast` for Stage 12,
`iqtree_ml_ufboot1000` for Stage 19. Two tests pin the command: one asserts
`--fast` and `--redo` are present for the fast build, the other asserts
`-B 1000 --bnni` are present and `--fast` is *absent* for the bootstrap build.
`[TEST]`

`tree_output_prefix` names outputs `<organelle>.<run_label>.iqtree_ml`, so the
cpDNA final tree is `cpDNA.primary.iqtree_ml.treefile`. A test asserts the exact
prefix for an mtDNA/primary input. `[TEST]`

### Resumability and the treefile guard

Like the other heavy stages, tree building reuses an existing `.treefile` unless
`--force` is passed. After the run, it re-checks that the treefile exists and is
non-empty, raising `PhylogeneticTreeError` otherwise — IQ-TREE exiting 0 is not
enough; the actual product must be there. `[CODE]`

## 11.4 Stage 14/20 — rendering

[`tree_visualization.py`](../tree_visualization.py) reads the tree summary,
parses each `.treefile` with `Bio.Phylo.read(path, "newick")`, and draws it with
matplotlib. It **does not alter the topology** — it only renders. Figure height
scales with the number of tips so a 275-tip tree is legible:

```python
def compute_tree_figure_size(sample_count: int) -> tuple[float, float]:
    width = 14.0
    height = max(6.0, min(80.0, sample_count * 0.16 + 2.0))
    return width, height
```

For 275 tips that is a tall 14×46-inch figure. Tests confirm the height grows
with sample count and that PNG/PDF/SVG are all written with a positive tip count.
`[TEST]` Stage 20 is the same code pointed at the Stage 19 trees, which is why
the final publication figures show bootstrap support and the Stage 14 figures do
not.

## 11.5 The Python concepts here

- **Boolean flag composition** (`fast = not A and not B`) to make two run modes
  from one function.
- **`shutil.which("iqtree") or shutil.which("iqtree2")`** to accept either
  executable name.
- **Biopython `Phylo`** for Newick parsing and drawing — the pipeline's one use
  of a real phylogenetics library.
- **Post-run output verification** (treefile exists and is non-empty) as a
  separate check from the subprocess return code.

## 11.6 The result, stated honestly

The final deliverable trees are the Stage 19 bootstrap trees, rendered in Stage
20, for 275 samples on 124,538 cpDNA sites and 44,930 mtDNA sites, under
`GTR+F+G4` with 1,000 ultrafast bootstrap replicates and BNNI. `[RESULT]`

What a tree here **is**: a maximum-likelihood estimate of relationships among the
sampled organelle sequences, with ultrafast-bootstrap support on internal
branches.

What it is **not**: a species tree, and not evidence of nuclear relationships.
Because each organelle is one largely linked cytoplasmic locus
([Chapter 6, §6.2](./06-organelle-biology.md)), even a fully bootstrap-supported
clade is strong evidence about *organelle-haplotype* history only. And the mtDNA
tree, built from only ~146 variable sites, has less site-level resolution than the
cpDNA tree despite having the same 275 tips. Its overall consensus missingness is
low (0.2534%), not a primary explanation for that difference. How to read support
values, cpDNA-vs-mtDNA discordance, and the limits of these trees is [Chapter
15](./15-phylogenetics-interpretation.md). `[BIO]`

## 11.7 Failure modes

- **Missing `iqtree`/`iqtree2`** → `PhylogeneticTreeError` before any run.
  `[CODE]`
- **Missing consensus FASTA** → raised by `read_tree_inputs`. `[CODE]`
- **IQ-TREE exits non-zero, or writes no treefile** → `PhylogeneticTreeError`
  naming the log. `[CODE]`
- **Missing plotting dependency (matplotlib/Biopython)** →
  `TreeVisualizationError` pointing at the tool audit. `[CODE]`
- **Scientific failure: reading the fast tree as final.** The Stage 12 tree has
  no bootstrap support; treating its topology as confident is a misread. Always
  confirm you are looking at Stage 19/20 for any claim. `[BIO]`

## 11.8 Exercises

1. **Trace.** With `--bootstrap-replicates 1000` and no `--full-search`, what is
   the value of `fast` in `main`, and which flags does `build_iqtree_command`
   append?
2. **Predict.** You run Stage 12 with `--full-search`. What `method` label is
   recorded, and is `--fast` in the command?
3. **Predict.** `compute_tree_figure_size(5)` versus `compute_tree_figure_size(275)`
   — give both heights and explain the `min(80.0, ...)` cap.
4. **Modify.** You want a third tree run with 5,000 UFBoot replicates in a new
   `results/21_.../` directory without touching Stage 19. What runner arguments
   do it, and which module code changes (if any)?
5. **Debug.** IQ-TREE returns 0 but the stage still raises
   `PhylogeneticTreeError`. What condition triggers that, and where would you
   look first?
6. **Interpret.** The cpDNA tree and mtDNA tree place the same population in
   different clades. Before calling this "discordance," name two data-quality
   explanations you must rule out using the Stage 11 and Stage 09 numbers.

Solutions in [Chapter 19](./19-solutions.md).

> Next: [Chapter 12 — PCA and Admixture-Style Clustering (Stages 15, 16, 18)](./12-pca-and-clustering.md)
