# Population Genetics

This step computes pairwise population Fst and per-population summary
statistics from filtered haploid cpDNA and mtDNA SNP alignments.
Only samples with resolved population codes are included in these
population-level summaries.

## Run

- Run label: `primary`
- Fst: Nei-style haploid SNP differentiation, averaged across informative sites
- Population summaries: sample count, haplotypes, diversity, nucleotide diversity, private variant sites

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Total samples in SNP alignment: 275
- Metadata-resolved populations: 34
- Pairwise comparisons: 561
- Pairwise Fst table: `dudleya_organelle_alignment_pipeline/results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv`
- Population summary table: `dudleya_organelle_alignment_pipeline/results/17_population_genetics/cpDNA.primary.population_genetics.population_summary.tsv`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Total samples in SNP alignment: 275
- Metadata-resolved populations: 34
- Pairwise comparisons: 561
- Pairwise Fst table: `dudleya_organelle_alignment_pipeline/results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv`
- Population summary table: `dudleya_organelle_alignment_pipeline/results/17_population_genetics/mtDNA.primary.population_genetics.population_summary.tsv`
