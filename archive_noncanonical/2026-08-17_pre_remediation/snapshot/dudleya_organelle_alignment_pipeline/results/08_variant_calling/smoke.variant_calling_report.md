# Haploid Variant Calling

This step calls raw haploid variants separately for cpDNA and mtDNA.
Filtering and consensus generation happen in later steps.

## Run

- Run label: `smoke`
- Samples called: 5
- Minimum mapping quality: 20
- Minimum base quality: 20
- Per-sample maximum pileup depth: 10000

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Raw variant records: 556
- Raw VCF: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/cpDNA.smoke.raw.vcf.gz`
- Index: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/cpDNA.smoke.raw.vcf.gz.tbi`
- Log: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/cpDNA.smoke.raw.bcftools.log`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Raw variant records: 40
- Raw VCF: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/mtDNA.smoke.raw.vcf.gz`
- Index: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/mtDNA.smoke.raw.vcf.gz.tbi`
- Log: `dudleya_organelle_alignment_pipeline/results/08_variant_calling/mtDNA.smoke.raw.bcftools.log`
