# Review Response Outputs

Derived outputs responding to the tree-rooting and PCA-legend requests.
The Stage 19 IQ-TREE searches were not rerun; those trees already include 1,000 ultrafast bootstrap replicates with BNNI correction.

## Publication Figures

- [Paper-style figure set](publication_figures/README.md)
- [Rooted cpDNA and mtDNA trees](publication_figures/figure_1_rooted_collapsed_trees.png)
- [cpDNA/mtDNA branch comparison](publication_figures/figure_2_cpdna_mtdna_branch_comparison.png)
- [Legended cpDNA and mtDNA PCA](publication_figures/figure_3_pca_requested_groups.png)

## Rooted Trees

- cpDNA rooted with ABAB + ABMU: [treefile](cpDNA.primary.rooted_ABAB_ABMU.iqtree_ml.treefile), [PNG](cpDNA.primary.rooted_ABAB_ABMU.iqtree_ml.png), [PDF](cpDNA.primary.rooted_ABAB_ABMU.iqtree_ml.pdf), [SVG](cpDNA.primary.rooted_ABAB_ABMU.iqtree_ml.svg)
- mtDNA rooted with ABAB + ABMU: [treefile](mtDNA.primary.rooted_ABAB_ABMU.iqtree_ml.treefile), [PNG](mtDNA.primary.rooted_ABAB_ABMU.iqtree_ml.png), [PDF](mtDNA.primary.rooted_ABAB_ABMU.iqtree_ml.pdf), [SVG](mtDNA.primary.rooted_ABAB_ABMU.iqtree_ml.svg)

Alternative ABAB-only and ABMU-only rooted versions are also included in `rooted_tree_summary.tsv`.

## Tree Comparison

- [cpDNA_mtDNA_tree_comparison.md](cpDNA_mtDNA_tree_comparison.md)
- [strongly_supported_shared_splits.tsv](strongly_supported_shared_splits.tsv)
- [strongly_supported_conflicting_splits.tsv](strongly_supported_conflicting_splits.tsv)

## PCA With Embedded Legend

- cpDNA PCA grouped by requested groups: [PNG](cpDNA.primary.pca.requested_groups.png), [PDF](cpDNA.primary.pca.requested_groups.pdf), [SVG](cpDNA.primary.pca.requested_groups.svg), [coordinates](cpDNA.primary.pca.requested_groups.coordinates.tsv)
- mtDNA PCA grouped by requested groups: [PNG](mtDNA.primary.pca.requested_groups.png), [PDF](mtDNA.primary.pca.requested_groups.pdf), [SVG](mtDNA.primary.pca.requested_groups.svg), [coordinates](mtDNA.primary.pca.requested_groups.coordinates.tsv)

Group mapping used here:

- `ABAB`: ABAB-prefixed samples / D. abramsii ssp. abramsii
- `ABBE`: ABBE-prefixed samples / D. abramsii ssp. bettinae
- `ABMU`: ABMU-prefixed samples / D. abramsii ssp. murina
- `DUSE`: samples annotated as D. setchellii
- `DUCY`: samples annotated as D. cymosa
- `Other / legacy IDs`: included samples without one of those labels
