# Variant Filtering

This step filters the raw haploid cpDNA and mtDNA variant calls.
Consensus FASTA generation, alignments, PCA, and trees happen in later steps.

## Run

- Run label: `primary`
- Variant type retained: biallelic SNPs
- Maximum missing genotype fraction: 0.2
- Minimum minor allele count: 2

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Samples: 278
- Raw records: 2531
- Filtered records: 2022
- Filtered VCF: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/09_variant_filtering/cpDNA.primary.filtered.vcf.gz`
- Index: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/09_variant_filtering/cpDNA.primary.filtered.vcf.gz.tbi`
- Log: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/09_variant_filtering/cpDNA.primary.filtered.bcftools.log`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Samples: 278
- Raw records: 192
- Filtered records: 146
- Filtered VCF: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/09_variant_filtering/mtDNA.primary.filtered.vcf.gz`
- Index: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/09_variant_filtering/mtDNA.primary.filtered.vcf.gz.tbi`
- Log: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/09_variant_filtering/mtDNA.primary.filtered.bcftools.log`
