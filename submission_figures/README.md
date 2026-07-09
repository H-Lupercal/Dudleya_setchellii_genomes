# Submission Figures

This folder contains extra visual summaries built from the completed
`full_pipeline_run/` outputs.

## Figures

| File | What it shows |
|---|---|
| [dudleya_organelle_submission_panel.png](dudleya_organelle_submission_panel.png) | Six-panel summary with cpDNA/mtDNA circular maps, PCA plots, and ADMIXTURE plots. |
| [cpDNA.pairwise_fst_heatmap.png](cpDNA.pairwise_fst_heatmap.png) | Pairwise cpDNA Fst among populations. |
| [mtDNA.pairwise_fst_heatmap.png](mtDNA.pairwise_fst_heatmap.png) | Pairwise mtDNA Fst among populations. |
| [cpDNA.population_colored_tree.png](cpDNA.population_colored_tree.png) | cpDNA maximum-likelihood tree with tip labels colored by population group. |
| [mtDNA.population_colored_tree.png](mtDNA.population_colored_tree.png) | mtDNA maximum-likelihood tree with tip labels colored by population group. |
| [submission_figures_summary.tsv](submission_figures_summary.tsv) | List of generated figure files. |

PDF and SVG versions are also written for the heatmaps and tree figures. The
multi-panel summary is written as PNG and PDF.

## Build

```bash
.tools/bioconda-env/bin/python3 submission_figures/build_submission_figures.py
```

The script uses the final rerun tables and figures under `full_pipeline_run/`
plus the circular genome maps under `genome_maps/`.
