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
- Mean CV error at best K: 0.08898600
- CV-error SD at best K: 0.01449154
- Replicates per K: 5
- Structure plot: `dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png`
- Best-K Q table: `dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/cpDNA.primary.bestK8.q.tsv`
- CV plot: `dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/cpDNA.primary.admixture_cv.png`

| K | Mean CV error | SD | Replicates | Best |
|---|---:|---:|---:|---|
| 1 | 0.48107400 | 0.00047595 | 5 | no |
| 2 | 0.28445000 | 0.00051648 | 5 | no |
| 3 | 0.23374800 | 0.00985176 | 5 | no |
| 4 | 0.18301800 | 0.00535419 | 5 | no |
| 5 | 0.16354200 | 0.00518844 | 5 | no |
| 6 | 0.12633000 | 0.00452729 | 5 | no |
| 7 | 0.09739600 | 0.00949484 | 5 | no |
| 8 | 0.08898600 | 0.01449154 | 5 | yes |

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Best K: 8
- Mean CV error at best K: 0.12644400
- CV-error SD at best K: 0.02207443
- Replicates per K: 5
- Structure plot: `dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png`
- Best-K Q table: `dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/mtDNA.primary.bestK8.q.tsv`
- CV plot: `dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/mtDNA.primary.admixture_cv.png`

| K | Mean CV error | SD | Replicates | Best |
|---|---:|---:|---:|---|
| 1 | 0.44209200 | 0.00153570 | 5 | no |
| 2 | 0.27713600 | 0.00151810 | 5 | no |
| 3 | 0.23636400 | 0.00825067 | 5 | no |
| 4 | 0.19291800 | 0.00365670 | 5 | no |
| 5 | 0.16777400 | 0.01279960 | 5 | no |
| 6 | 0.14287400 | 0.01610858 | 5 | no |
| 7 | 0.12894200 | 0.01863814 | 5 | no |
| 8 | 0.12644400 | 0.02207443 | 5 | yes |
