# Admixture-Style Clustering

This step runs ADMIXTURE on cpDNA and mtDNA SNP alignments separately.
Because ADMIXTURE is a diploid-oriented tool, haploid organelle calls
are encoded as pseudo-diploid homozygous genotypes. These plots should
be interpreted as organelle haplotype clustering, not nuclear admixture.

## Run

- Run label: `primary`
- K selection: lowest cross-validation error among tested K values
- Haploid encoding: each called base is duplicated; missing calls are `0 0`

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Best K: 8
- Best CV error: 0.07477000
- Structure plot: `dudleya_organelle_alignment_pipeline/results/16_admixture/cpDNA.primary.bestK8.structure.png`
- Best-K Q table: `dudleya_organelle_alignment_pipeline/results/16_admixture/cpDNA.primary.bestK8.q.tsv`
- CV plot: `dudleya_organelle_alignment_pipeline/results/16_admixture/cpDNA.primary.admixture_cv.png`

| K | CV error | Best |
|---|---:|---|
| 1 | 0.48074000 | no |
| 2 | 0.28426000 | no |
| 3 | 0.22297000 | no |
| 4 | 0.18187000 | no |
| 5 | 0.16561000 | no |
| 6 | 0.12722000 | no |
| 7 | 0.08946000 | no |
| 8 | 0.07477000 | yes |

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Best K: 7
- Best CV error: 0.11670000
- Structure plot: `dudleya_organelle_alignment_pipeline/results/16_admixture/mtDNA.primary.bestK7.structure.png`
- Best-K Q table: `dudleya_organelle_alignment_pipeline/results/16_admixture/mtDNA.primary.bestK7.q.tsv`
- CV plot: `dudleya_organelle_alignment_pipeline/results/16_admixture/mtDNA.primary.admixture_cv.png`

| K | CV error | Best |
|---|---:|---|
| 1 | 0.44062000 | no |
| 2 | 0.27600000 | no |
| 3 | 0.23043000 | no |
| 4 | 0.19754000 | no |
| 5 | 0.15383000 | no |
| 6 | 0.15097000 | no |
| 7 | 0.11670000 | yes |
| 8 | 0.14368000 | no |
