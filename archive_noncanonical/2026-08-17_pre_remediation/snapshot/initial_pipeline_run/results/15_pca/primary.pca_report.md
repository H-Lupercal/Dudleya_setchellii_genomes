# PCA Visualization

This step computes cpDNA and mtDNA PCA from the filtered haploid
SNP-only alignments. Missing SNP states are mean-imputed per retained
site before PCA, and plots are colored by available species/population
metadata.

## Run

- Run label: `primary`
- Input: Step 9 SNP-only FASTA alignments
- Output formats: coordinates TSV, variance TSV, PNG, PDF, SVG

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Samples: 275
- SNP alignment sites: 2015
- Retained polymorphic sites: 2015
- PC1 variance: 36.62%
- PC2 variance: 14.65%
- Coordinates: `dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.coordinates.tsv`
- PNG: `dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.png`
- PDF: `dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.pdf`
- SVG: `dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.svg`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Samples: 275
- SNP alignment sites: 146
- Retained polymorphic sites: 146
- PC1 variance: 34.48%
- PC2 variance: 14.06%
- Coordinates: `dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.coordinates.tsv`
- PNG: `dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.png`
- PDF: `dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.pdf`
- SVG: `dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.svg`
