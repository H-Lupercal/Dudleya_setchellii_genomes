# Initial Pipeline Run

This folder keeps the final-output snapshot from the first organelle population
genomics run.

The files under `results/` are copied from the original pipeline results
directory so they can be browsed directly on GitHub. The original location is
also kept as `source_results_symlink`:

```text
initial_pipeline_run/source_results_symlink -> ../dudleya_organelle_alignment_pipeline/results
```

## Final Results

| File | What it is |
|---|---|
| [results/organelle_population_report.md](results/organelle_population_report.md) | Initial run population-genomics report. |
| [results/final_deliverables_manifest.tsv](results/final_deliverables_manifest.tsv) | Initial run final deliverables manifest. |
| [results/10_snp_alignment/cpDNA.primary.snp_alignment.fa](results/10_snp_alignment/cpDNA.primary.snp_alignment.fa) | Initial cpDNA SNP alignment. |
| [results/10_snp_alignment/mtDNA.primary.snp_alignment.fa](results/10_snp_alignment/mtDNA.primary.snp_alignment.fa) | Initial mtDNA SNP alignment. |
| [results/11_callable_consensus/cpDNA.primary.callable_consensus.fa](results/11_callable_consensus/cpDNA.primary.callable_consensus.fa) | Initial cpDNA callable consensus alignment. |
| [results/11_callable_consensus/mtDNA.primary.callable_consensus.fa](results/11_callable_consensus/mtDNA.primary.callable_consensus.fa) | Initial mtDNA callable consensus alignment. |
| [results/15_pca/cpDNA.primary.pca.png](results/15_pca/cpDNA.primary.pca.png) | Initial cpDNA PCA figure. |
| [results/15_pca/mtDNA.primary.pca.png](results/15_pca/mtDNA.primary.pca.png) | Initial mtDNA PCA figure. |
| [results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png](results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png) | Initial rendered cpDNA ML tree figure. |
| [results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png](results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png) | Initial rendered mtDNA ML tree figure. |

Use [../full_pipeline_run/](../full_pipeline_run/) for the main 16-thread rerun.
