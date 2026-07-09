# Second Full Pipeline Run Results

This folder contains the second complete rerun of the Dudleya cpDNA/mtDNA
organelle population-genomics pipeline.

Run facts:

- Run directory: `full_pipeline_run/`
- Started: `2026-07-08T09:24:13-07:00`
- Finished: `2026-07-08T18:46:47-07:00`
- CPU threads requested: `16`
- Result size: `34G`
- Final downstream sample count: `278`
- Stage ledger: [`logs/stage_status.tsv`](logs/stage_status.tsv)
- Run metadata: [`run_metadata.txt`](run_metadata.txt)
- Full narrative report: [`FULL_PIPELINE_RUN_REPORT.md`](FULL_PIPELINE_RUN_REPORT.md)

## Folder Layout

| Path | Contents |
|---|---|
| [`results/`](results/) | Organized per-stage outputs from `00_manifest` through `20_bootstrap_tree_visualization`. |
| [`logs/`](logs/) | One log per pipeline stage plus the stage-completion ledger. |
| [`run_full_pipeline.sh`](run_full_pipeline.sh) | Start-to-finish runner used for this rerun. |
| [`resume_from_stage_07.sh`](resume_from_stage_07.sh) | Resume helper used after the sample-count guard was updated. |
| [`resume_from_stage_15.sh`](resume_from_stage_15.sh) | Resume helper used after selecting the bioconda Python environment. |
| [`resume_from_stage_16.sh`](resume_from_stage_16.sh) | Resume helper used after ADMIXTURE all-missing-sample handling was added. |

## Final Reports And Ledgers

| File | What it is |
|---|---|
| [`FULL_PIPELINE_RUN_REPORT.md`](FULL_PIPELINE_RUN_REPORT.md) | Human-readable report covering what ran, what changed, final outputs, and verification. |
| [`logs/stage_status.tsv`](logs/stage_status.tsv) | Completion ledger for all stages `00` through `20`. |
| [`results/13_tool_audit/primary.tool_audit_report.md`](results/13_tool_audit/primary.tool_audit_report.md) | Tool audit from the rerun environment. |
| [`results/18_admixture_replicates/primary.admixture_report.md`](results/18_admixture_replicates/primary.admixture_report.md) | Five-replicate ADMIXTURE-style clustering report. |
| [`results/19_bootstrap_phylogenetic_tree/primary.phylogenetic_tree_report.md`](results/19_bootstrap_phylogenetic_tree/primary.phylogenetic_tree_report.md) | 1,000-UFBoot IQ-TREE report. |
| [`results/20_bootstrap_tree_visualization/primary.tree_visualization_report.md`](results/20_bootstrap_tree_visualization/primary.tree_visualization_report.md) | Final bootstrap tree rendering report. |

## Final Alignments

| File | What it is |
|---|---|
| [`results/11_callable_consensus/cpDNA.primary.callable_consensus.fa`](results/11_callable_consensus/cpDNA.primary.callable_consensus.fa) | cpDNA callable-site consensus alignment: 278 samples x 124,538 sites. |
| [`results/11_callable_consensus/mtDNA.primary.callable_consensus.fa`](results/11_callable_consensus/mtDNA.primary.callable_consensus.fa) | mtDNA callable-site consensus alignment: 278 samples x 44,930 sites. |
| [`results/10_snp_alignment/cpDNA.primary.snp_alignment.fa`](results/10_snp_alignment/cpDNA.primary.snp_alignment.fa) | cpDNA filtered haploid SNP alignment: 278 samples x 2,022 SNP sites. |
| [`results/10_snp_alignment/mtDNA.primary.snp_alignment.fa`](results/10_snp_alignment/mtDNA.primary.snp_alignment.fa) | mtDNA filtered haploid SNP alignment: 278 samples x 146 SNP sites. |
| [`results/10_snp_alignment/primary.snp_alignment_summary.tsv`](results/10_snp_alignment/primary.snp_alignment_summary.tsv) | Machine-readable SNP alignment summary. |
| [`results/11_callable_consensus/primary.callable_consensus_summary.tsv`](results/11_callable_consensus/primary.callable_consensus_summary.tsv) | Machine-readable callable-consensus summary. |

## PCA

| File | What it is |
|---|---|
| [`results/15_pca/cpDNA.primary.pca.png`](results/15_pca/cpDNA.primary.pca.png) | cpDNA PCA plot; PC1 = 37.04 percent, PC2 = 14.45 percent. |
| [`results/15_pca/mtDNA.primary.pca.png`](results/15_pca/mtDNA.primary.pca.png) | mtDNA PCA plot; PC1 = 34.43 percent, PC2 = 14.03 percent. |
| [`results/15_pca/cpDNA.primary.pca.coordinates.tsv`](results/15_pca/cpDNA.primary.pca.coordinates.tsv) | cpDNA PCA coordinates and metadata. |
| [`results/15_pca/mtDNA.primary.pca.coordinates.tsv`](results/15_pca/mtDNA.primary.pca.coordinates.tsv) | mtDNA PCA coordinates and metadata. |
| [`results/15_pca/primary.pca_summary.tsv`](results/15_pca/primary.pca_summary.tsv) | PCA summary table. |

## Phylogenetic Trees

| File | What it is |
|---|---|
| [`results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile`](results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile) | cpDNA IQ-TREE maximum-likelihood tree with 1,000 ultrafast bootstraps and BNNI. |
| [`results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile`](results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile) | mtDNA IQ-TREE maximum-likelihood tree with 1,000 ultrafast bootstraps and BNNI. |
| [`results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png`](results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png) | Rendered cpDNA bootstrap tree figure. |
| [`results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png`](results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png) | Rendered mtDNA bootstrap tree figure. |
| [`results/19_bootstrap_phylogenetic_tree/primary.phylogenetic_tree_summary.tsv`](results/19_bootstrap_phylogenetic_tree/primary.phylogenetic_tree_summary.tsv) | Bootstrap-tree summary table. |
| [`results/20_bootstrap_tree_visualization/primary.tree_visualization_summary.tsv`](results/20_bootstrap_tree_visualization/primary.tree_visualization_summary.tsv) | Tree-figure summary table. |

## ADMIXTURE And Population Statistics

| File | What it is |
|---|---|
| [`results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png`](results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png) | cpDNA five-replicate ADMIXTURE-style structure plot; best K = 8 by mean CV error. |
| [`results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png`](results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png) | mtDNA five-replicate ADMIXTURE-style structure plot; best K = 8 by mean CV error. |
| [`results/18_admixture_replicates/primary.admixture_summary.tsv`](results/18_admixture_replicates/primary.admixture_summary.tsv) | K=1..8, five-replicate CV-error summary. |
| [`results/18_admixture_replicates/mtDNA.primary.pseudo_diploid.excluded_samples.tsv`](results/18_admixture_replicates/mtDNA.primary.pseudo_diploid.excluded_samples.tsv) | One mtDNA ADMIXTURE-only exclusion for an all-missing SNP genotype sample. |
| [`results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv`](results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv) | cpDNA pairwise Fst table: 595 comparisons across 35 populations. |
| [`results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv`](results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv) | mtDNA pairwise Fst table: 595 comparisons across 35 populations. |
| [`results/17_population_genetics/cpDNA.primary.population_genetics.population_summary.tsv`](results/17_population_genetics/cpDNA.primary.population_genetics.population_summary.tsv) | cpDNA population summary table. |
| [`results/17_population_genetics/mtDNA.primary.population_genetics.population_summary.tsv`](results/17_population_genetics/mtDNA.primary.population_genetics.population_summary.tsv) | mtDNA population summary table. |
