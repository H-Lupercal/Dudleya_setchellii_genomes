# Second Full Pipeline Run Report

This report documents the second complete rerun of the Dudleya cpDNA/mtDNA
organelle population-genomics pipeline in [`full_pipeline_run/`](.).

## Executive Summary

The pipeline completed all stages from `00_manifest` through
`20_bootstrap_tree_visualization` using 16 CPU threads where pipeline tools
accepted a thread count. Outputs are organized under [`results/`](results/) and
stage logs are under [`logs/`](logs/).

Final run facts:

- Started: `2026-07-08T09:24:13-07:00`
- Finished: `2026-07-08T18:46:47-07:00`
- Result size: `34G`
- Final downstream sample count: `278`
- Stage ledger: [`logs/stage_status.tsv`](logs/stage_status.tsv)
- Run metadata: [`run_metadata.txt`](run_metadata.txt)

## What Happened During The Run

1. The full pipeline was launched with `16` requested CPU threads.
2. Stages `00` through `06` completed, producing manifests, reference checks,
   pilot alignments, mask definitions, and all-sample BAM/QC outputs.
3. Stage `07` initially stopped because the old guard expected `275`
   downstream samples while the rerun produced `278`. The guard was changed to
   accept the actual included sample set, and the run resumed from Stage `07`.
4. Stages `08` through `14` completed, including haploid variant calling,
   variant filtering, SNP alignments, callable consensus alignments, initial ML
   trees, tool audit, and initial tree figures.
5. Stage `15` initially used system Python, which did not have `scikit-learn`.
   The runner was updated to use the local bioconda Python at
   `.tools/bioconda-env/bin/python3`, then PCA completed.
6. Stage `16` initially hit an ADMIXTURE input error because one mtDNA sample
   had all SNP genotypes missing. The ADMIXTURE stage was updated to exclude
   all-missing samples from ADMIXTURE inputs only and record the exclusion.
7. Stages `16` through `20` then completed: single-run ADMIXTURE, population
   genetics, five-replicate ADMIXTURE, 1,000-UFBoot ML trees, and final tree
   visualizations.

## Final Deliverables

### Alignments

| Output | Result |
|---|---|
| [`results/11_callable_consensus/cpDNA.primary.callable_consensus.fa`](results/11_callable_consensus/cpDNA.primary.callable_consensus.fa) | 278 samples x 124,538 callable cpDNA sites. |
| [`results/11_callable_consensus/mtDNA.primary.callable_consensus.fa`](results/11_callable_consensus/mtDNA.primary.callable_consensus.fa) | 278 samples x 44,930 callable mtDNA sites. |
| [`results/10_snp_alignment/cpDNA.primary.snp_alignment.fa`](results/10_snp_alignment/cpDNA.primary.snp_alignment.fa) | 278 samples x 2,022 cpDNA SNP sites. |
| [`results/10_snp_alignment/mtDNA.primary.snp_alignment.fa`](results/10_snp_alignment/mtDNA.primary.snp_alignment.fa) | 278 samples x 146 mtDNA SNP sites. |

### PCA

| Output | Result |
|---|---|
| [`results/15_pca/cpDNA.primary.pca.png`](results/15_pca/cpDNA.primary.pca.png) | cpDNA PCA; PC1 = 37.04 percent, PC2 = 14.45 percent. |
| [`results/15_pca/mtDNA.primary.pca.png`](results/15_pca/mtDNA.primary.pca.png) | mtDNA PCA; PC1 = 34.43 percent, PC2 = 14.03 percent. |

### Phylogenetic Trees

| Output | Result |
|---|---|
| [`results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile`](results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile) | cpDNA IQ-TREE ML tree, GTR+F+G4, 1,000 ultrafast bootstraps with BNNI. |
| [`results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile`](results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile) | mtDNA IQ-TREE ML tree, GTR+F+G4, 1,000 ultrafast bootstraps with BNNI. |
| [`results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png`](results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png) | Rendered cpDNA bootstrap tree, 278 tips. |
| [`results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png`](results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png) | Rendered mtDNA bootstrap tree, 278 tips. |

### ADMIXTURE And Population Genetics

| Output | Result |
|---|---|
| [`results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png`](results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png) | cpDNA best K = 8; five-replicate mean CV error = 0.08512200. |
| [`results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png`](results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png) | mtDNA best K = 8; five-replicate mean CV error = 0.11156800. |
| [`results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv`](results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv) | cpDNA pairwise Fst: 595 comparisons across 35 populations. |
| [`results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv`](results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv) | mtDNA pairwise Fst: 595 comparisons across 35 populations. |

## ADMIXTURE-Only Exclusion

ADMIXTURE rejects individuals with all genotypes missing. One mtDNA sample was
therefore excluded from ADMIXTURE input only:

| Sample | Organelle | Reason |
|---|---|---|
| `CY_RED_LP_202_Du-561` | mtDNA | `all_snp_genotypes_missing` |

The exclusion file is
[`results/18_admixture_replicates/mtDNA.primary.pseudo_diploid.excluded_samples.tsv`](results/18_admixture_replicates/mtDNA.primary.pseudo_diploid.excluded_samples.tsv).
Other pipeline stages retain the full 278-sample downstream set.

## Verification

The final checks were:

- All stages `00` through `20` appear in
  [`logs/stage_status.tsv`](logs/stage_status.tsv).
- Final tree visualization rendered cpDNA and mtDNA trees with 278 tips each.
- Unit tests passed: `70 tests in 1.149s, OK`.

## Notes And Caveats

- ADMIXTURE is used as an organelle haplotype-clustering visualization. Haploid
  organelle calls are encoded as pseudo-diploid homozygous genotypes.
- mtDNA population-genetic interpretation uses the high-confidence unique track
  rather than repeat-rich mtDNA regions.
- The final trees are maximum-likelihood IQ-TREE outputs with ultrafast
  bootstrap support, not neighbor-joining trees.
