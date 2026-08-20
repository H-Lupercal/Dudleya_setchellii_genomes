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
- Total samples in SNP alignment: 278
- Metadata-resolved populations: 35
- Pairwise comparisons: 595
- Pairwise Fst table: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv`
- Population summary table: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/17_population_genetics/cpDNA.primary.population_genetics.population_summary.tsv`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Total samples in SNP alignment: 278
- Metadata-resolved populations: 35
- Pairwise comparisons: 595
- Pairwise Fst table: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv`
- Population summary table: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/17_population_genetics/mtDNA.primary.population_genetics.population_summary.tsv`
