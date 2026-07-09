# Variant Filtering

This step filters the raw haploid cpDNA and mtDNA variant calls from Step 7.
Consensus FASTA generation, alignments, PCA, and trees happen in later steps.

## Run

- Run label: `primary`
- Variant type retained: biallelic SNPs
- Maximum missing genotype fraction: 0.2
- Minimum minor allele count: 2

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Samples: 275
- Raw records: 2475
- Filtered records: 2015
- Filtered VCF: `dudleya_organelle_alignment_pipeline/results/09_variant_filtering/cpDNA.primary.filtered.vcf.gz`
- Index: `dudleya_organelle_alignment_pipeline/results/09_variant_filtering/cpDNA.primary.filtered.vcf.gz.tbi`
- Log: `dudleya_organelle_alignment_pipeline/results/09_variant_filtering/cpDNA.primary.filtered.bcftools.log`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Samples: 275
- Raw records: 190
- Filtered records: 146
- Filtered VCF: `dudleya_organelle_alignment_pipeline/results/09_variant_filtering/mtDNA.primary.filtered.vcf.gz`
- Index: `dudleya_organelle_alignment_pipeline/results/09_variant_filtering/mtDNA.primary.filtered.vcf.gz.tbi`
- Log: `dudleya_organelle_alignment_pipeline/results/09_variant_filtering/mtDNA.primary.filtered.bcftools.log`
