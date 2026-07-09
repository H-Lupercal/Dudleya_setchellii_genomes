# Admixture-Style Clustering

This step runs ADMIXTURE on cpDNA and mtDNA SNP alignments separately.
Because ADMIXTURE is a diploid-oriented tool, haploid organelle calls
are encoded as pseudo-diploid homozygous genotypes. These plots should
be interpreted as organelle haplotype clustering, not nuclear admixture.

## Run

- Run label: `primary`
- K selection: lowest mean cross-validation error among tested K values
- Haploid encoding: each called base is duplicated; missing calls are `0 0`

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Best K: 8
- Mean CV error at best K: 0.08512200
- CV-error SD at best K: 0.01378138
- Replicates per K: 5
- Structure plot: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png`
- Best-K Q table: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/18_admixture_replicates/cpDNA.primary.bestK8.q.tsv`
- CV plot: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/18_admixture_replicates/cpDNA.primary.admixture_cv.png`

| K | Mean CV error | SD | Replicates | Best |
|---|---:|---:|---:|---|
| 1 | 0.48286800 | 0.00038868 | 5 | no |
| 2 | 0.28462400 | 0.00049933 | 5 | no |
| 3 | 0.23357800 | 0.01321700 | 5 | no |
| 4 | 0.18761400 | 0.00389663 | 5 | no |
| 5 | 0.16149600 | 0.01307840 | 5 | no |
| 6 | 0.13668200 | 0.01372450 | 5 | no |
| 7 | 0.10698200 | 0.00985610 | 5 | no |
| 8 | 0.08512200 | 0.01378138 | 5 | yes |

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Best K: 8
- Mean CV error at best K: 0.11156800
- CV-error SD at best K: 0.00604363
- Replicates per K: 5
- Structure plot: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png`
- Best-K Q table: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/18_admixture_replicates/mtDNA.primary.bestK8.q.tsv`
- CV plot: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/18_admixture_replicates/mtDNA.primary.admixture_cv.png`

| K | Mean CV error | SD | Replicates | Best |
|---|---:|---:|---:|---|
| 1 | 0.44336400 | 0.00332705 | 5 | no |
| 2 | 0.27600800 | 0.00195051 | 5 | no |
| 3 | 0.23183200 | 0.00512357 | 5 | no |
| 4 | 0.19303600 | 0.00441710 | 5 | no |
| 5 | 0.16663400 | 0.00675772 | 5 | no |
| 6 | 0.14348800 | 0.01057158 | 5 | no |
| 7 | 0.13310200 | 0.01247851 | 5 | no |
| 8 | 0.11156800 | 0.00604363 | 5 | yes |
