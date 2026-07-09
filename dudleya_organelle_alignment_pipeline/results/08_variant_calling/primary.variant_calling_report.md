# Haploid Variant Calling

This step calls raw haploid variants separately for cpDNA and mtDNA.
Filtering and consensus generation happen in later steps.

## Run

- Run label: `primary`
- Samples called: 275
- Minimum mapping quality: 20
- Minimum base quality: 20
- Per-sample maximum pileup depth: 10000

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Raw variant records: 2475
- Raw VCF: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/cpDNA.primary.raw.vcf.gz`
- Index: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/cpDNA.primary.raw.vcf.gz.tbi`
- Log: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/cpDNA.primary.raw.bcftools.log`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Raw variant records: 190
- Raw VCF: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/mtDNA.primary.raw.vcf.gz`
- Index: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/mtDNA.primary.raw.vcf.gz.tbi`
- Log: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/mtDNA.primary.raw.bcftools.log`
