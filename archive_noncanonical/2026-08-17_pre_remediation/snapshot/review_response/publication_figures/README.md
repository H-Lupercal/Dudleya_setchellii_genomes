# Publication Figures

Static figures prepared from the rooted Stage 19 IQ-TREE results and Stage 15 PCA coordinates.
PDF and SVG files are the primary vector outputs; PNG files are 600 dpi previews.

| Figure | Purpose | Files |
|---|---|---|
| Figure 1 | Side-by-side cpDNA and mtDNA rooted trees. Colored points identify sample groups, black nodes mark UFBoot support of at least 95, and monophyletic population clades are collapsed where possible. | [PDF](figure_1_rooted_collapsed_trees.pdf), [SVG](figure_1_rooted_collapsed_trees.svg), [PNG](figure_1_rooted_collapsed_trees.png) |
| Figure 2 | Short comparison of strongly supported shared and incompatible cpDNA/mtDNA branches. | [PDF](figure_2_cpdna_mtdna_branch_comparison.pdf), [SVG](figure_2_cpdna_mtdna_branch_comparison.svg), [PNG](figure_2_cpdna_mtdna_branch_comparison.png) |
| Figure 3 | cpDNA and mtDNA PCA with a shared group legend, population centroids, and selected outlier labels. | [PDF](figure_3_pca_requested_groups.pdf), [SVG](figure_3_pca_requested_groups.svg), [PNG](figure_3_pca_requested_groups.png) |
| Supplementary Figure 1 | Complete cpDNA tree with all 278 sample labels. | [PDF](supplementary_figure_1_cpDNA_full_tree.pdf), [SVG](supplementary_figure_1_cpDNA_full_tree.svg), [PNG](supplementary_figure_1_cpDNA_full_tree.png) |
| Supplementary Figure 2 | Complete mtDNA tree with all 278 sample labels. | [PDF](supplementary_figure_2_mtDNA_full_tree.pdf), [SVG](supplementary_figure_2_mtDNA_full_tree.svg), [PNG](supplementary_figure_2_mtDNA_full_tree.png) |

The original Newick trees and PCA coordinate tables were not modified. Figure dimensions and formats are recorded in [figure_manifest.tsv](figure_manifest.tsv), and [collapsed_clades.tsv](collapsed_clades.tsv) lists every population clade collapsed in Figure 1.

Regenerate the figures from the repository root with:

```bash
.tools/bioconda-env/bin/Rscript review_response/build_publication_figures.R
```

