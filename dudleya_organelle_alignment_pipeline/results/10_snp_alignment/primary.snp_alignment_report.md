# SNP Alignment

This step converts filtered haploid cpDNA and mtDNA SNP VCFs into
SNP-only FASTA alignments. These alignments are intended for quick
tree-building and matrix-based analyses; full reference-length
consensus FASTAs can be generated in a later step if needed.

## Run

- Run label: `primary`
- Haploid genotype encoding: `0` uses REF, `1` uses ALT, missing uses `N`.

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Samples: 275
- Filtered records: 2015
- Alignment sites: 2015
- Missing alignment bases: 252
- FASTA: `dudleya_organelle_alignment_pipeline/results/10_snp_alignment/cpDNA.primary.snp_alignment.fa`
- Site table: `dudleya_organelle_alignment_pipeline/results/10_snp_alignment/cpDNA.primary.snp_sites.tsv`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Samples: 275
- Filtered records: 146
- Alignment sites: 146
- Missing alignment bases: 349
- FASTA: `dudleya_organelle_alignment_pipeline/results/10_snp_alignment/mtDNA.primary.snp_alignment.fa`
- Site table: `dudleya_organelle_alignment_pipeline/results/10_snp_alignment/mtDNA.primary.snp_sites.tsv`
