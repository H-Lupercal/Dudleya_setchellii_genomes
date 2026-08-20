# Response To Tree And PCA Questions

## 1. Rooting The Trees

Yes. The Stage 19 cpDNA and mtDNA maximum-likelihood trees were rerooted using
the ABAB and ABMU samples as the outgroup set. I also included ABAB-only and
ABMU-only rooted versions in case one outgroup choice is preferred during
review.

Primary rooted files:

- cpDNA ABAB+ABMU rooted tree: [treefile](cpDNA.primary.rooted_ABAB_ABMU.iqtree_ml.treefile), [PNG](cpDNA.primary.rooted_ABAB_ABMU.iqtree_ml.png), [PDF](cpDNA.primary.rooted_ABAB_ABMU.iqtree_ml.pdf), [SVG](cpDNA.primary.rooted_ABAB_ABMU.iqtree_ml.svg)
- mtDNA ABAB+ABMU rooted tree: [treefile](mtDNA.primary.rooted_ABAB_ABMU.iqtree_ml.treefile), [PNG](mtDNA.primary.rooted_ABAB_ABMU.iqtree_ml.png), [PDF](mtDNA.primary.rooted_ABAB_ABMU.iqtree_ml.pdf), [SVG](mtDNA.primary.rooted_ABAB_ABMU.iqtree_ml.svg)

The full list of rooted outputs and outgroup tips is in
[rooted_tree_summary.tsv](rooted_tree_summary.tsv).

## 2. Bootstrap Replicates And cpDNA/mtDNA Tree Comparison

The existing Stage 19 trees already include 1,000 ultrafast bootstrap replicates
with BNNI correction, so I did not rerun IQ-TREE. The treefiles with support
values are:

- [cpDNA.primary.iqtree_ml.treefile](../full_pipeline_run/results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile)
- [mtDNA.primary.iqtree_ml.treefile](../full_pipeline_run/results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile)

I added a short comparison focused on branches with UFBoot support >= 95:

- [cpDNA_mtDNA_tree_comparison.md](cpDNA_mtDNA_tree_comparison.md)
- [strongly_supported_shared_splits.tsv](strongly_supported_shared_splits.tsv)
- [strongly_supported_conflicting_splits.tsv](strongly_supported_conflicting_splits.tsv)

Short version: several strongly supported groups are shared between cpDNA and
mtDNA, including ABAB/ABMU outgroup structure, several DUCY groups, and several
DUSE groups. The clearest discrepancies are strongly supported DUSE placements
that differ between cpDNA and mtDNA, plus one DUCY-centered split where the same
general sample set is involved but one included CY_HICN sample differs between
the two trees.

## 3. PCA Legend And Group Colors

Yes. I regenerated the PCA plots with the legend embedded in the figure and
colored samples by the requested broad groups:

- `DUSE`: samples annotated as *D. setchellii*
- `DUCY`: samples annotated as *D. cymosa*
- `ABAB`: ABAB-prefixed samples / *D. abramsii* ssp. *abramsii*
- `ABBE`: ABBE-prefixed samples / *D. abramsii* ssp. *bettinae*
- `ABMU`: ABMU-prefixed samples / *D. abramsii* ssp. *murina*
- `Other / legacy IDs`: included samples without one of those labels

Updated PCA files:

- cpDNA PCA with embedded legend: [PNG](cpDNA.primary.pca.requested_groups.png), [PDF](cpDNA.primary.pca.requested_groups.pdf), [SVG](cpDNA.primary.pca.requested_groups.svg), [coordinates](cpDNA.primary.pca.requested_groups.coordinates.tsv)
- mtDNA PCA with embedded legend: [PNG](mtDNA.primary.pca.requested_groups.png), [PDF](mtDNA.primary.pca.requested_groups.pdf), [SVG](mtDNA.primary.pca.requested_groups.svg), [coordinates](mtDNA.primary.pca.requested_groups.coordinates.tsv)
