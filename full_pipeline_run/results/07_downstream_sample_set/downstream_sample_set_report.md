# Downstream Sample Set

This step defines the sample set for downstream haploid variant calling,
consensus FASTA generation, cpDNA/mtDNA all-sample alignments, PCA,
phylogenetic trees, Fst, and structure/admixture-style clustering.

## Summary

- Included samples: 278
- Excluded samples: 2

## Exclusions By Stage

- step0_manifest: 2

## Outputs

- `included_samples.tsv`: samples to use in primary downstream analyses.
- `excluded_samples.tsv`: samples excluded before downstream analyses.

The included sample set should be used by variant calling and
all later population-genetic outputs.
