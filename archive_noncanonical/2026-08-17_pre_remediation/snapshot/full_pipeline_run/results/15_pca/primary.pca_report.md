# PCA Visualization

This step computes cpDNA and mtDNA PCA from the filtered haploid
SNP-only alignments. Missing SNP states are mean-imputed per retained
site before PCA, and plots are colored by available species/population
metadata.

## Run

- Run label: `primary`
- Input: SNP-only FASTA alignments
- Output formats: coordinates TSV, variance TSV, PNG, PDF, SVG

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Samples: 278
- SNP alignment sites: 2022
- Retained polymorphic sites: 2022
- PC1 variance: 37.04%
- PC2 variance: 14.45%
- Coordinates: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/15_pca/cpDNA.primary.pca.coordinates.tsv`
- PNG: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/15_pca/cpDNA.primary.pca.png`
- PDF: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/15_pca/cpDNA.primary.pca.pdf`
- SVG: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/15_pca/cpDNA.primary.pca.svg`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Samples: 278
- SNP alignment sites: 146
- Retained polymorphic sites: 146
- PC1 variance: 34.43%
- PC2 variance: 14.03%
- Coordinates: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/15_pca/mtDNA.primary.pca.coordinates.tsv`
- PNG: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/15_pca/mtDNA.primary.pca.png`
- PDF: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/15_pca/mtDNA.primary.pca.pdf`
- SVG: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/15_pca/mtDNA.primary.pca.svg`
