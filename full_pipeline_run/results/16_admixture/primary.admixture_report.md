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
- Mean CV error at best K: 0.07502000
- CV-error SD at best K: 0.00000000
- Replicates per K: 1
- Structure plot: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/16_admixture/cpDNA.primary.bestK8.structure.png`
- Best-K Q table: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/16_admixture/cpDNA.primary.bestK8.q.tsv`
- CV plot: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/16_admixture/cpDNA.primary.admixture_cv.png`

| K | Mean CV error | SD | Replicates | Best |
|---|---:|---:|---:|---|
| 1 | 0.48327000 | 0.00000000 | 1 | no |
| 2 | 0.28453000 | 0.00000000 | 1 | no |
| 3 | 0.22632000 | 0.00000000 | 1 | no |
| 4 | 0.18503000 | 0.00000000 | 1 | no |
| 5 | 0.14261000 | 0.00000000 | 1 | no |
| 6 | 0.13212000 | 0.00000000 | 1 | no |
| 7 | 0.10782000 | 0.00000000 | 1 | no |
| 8 | 0.07502000 | 0.00000000 | 1 | yes |

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Best K: 8
- Mean CV error at best K: 0.11165000
- CV-error SD at best K: 0.00000000
- Replicates per K: 1
- Structure plot: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/16_admixture/mtDNA.primary.bestK8.structure.png`
- Best-K Q table: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/16_admixture/mtDNA.primary.bestK8.q.tsv`
- CV plot: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/16_admixture/mtDNA.primary.admixture_cv.png`

| K | Mean CV error | SD | Replicates | Best |
|---|---:|---:|---:|---|
| 1 | 0.44182000 | 0.00000000 | 1 | no |
| 2 | 0.27529000 | 0.00000000 | 1 | no |
| 3 | 0.23255000 | 0.00000000 | 1 | no |
| 4 | 0.19209000 | 0.00000000 | 1 | no |
| 5 | 0.16558000 | 0.00000000 | 1 | no |
| 6 | 0.12536000 | 0.00000000 | 1 | no |
| 7 | 0.13750000 | 0.00000000 | 1 | no |
| 8 | 0.11165000 | 0.00000000 | 1 | yes |
