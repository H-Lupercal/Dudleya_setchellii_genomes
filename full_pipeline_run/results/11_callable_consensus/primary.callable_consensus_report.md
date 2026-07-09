# Callable-Site Consensus Alignment

This step builds full callable-site FASTA alignments for cpDNA and mtDNA.
Each alignment follows the population-genetic BED track, starts
from the annotated organelle reference, overlays the filtered haploid
SNP genotypes, masks the raw variant sites that failed filtering, and
uses the depth files to write `N` at bases below the minimum depth.

## Run

- Run label: `primary`
- Minimum depth for a non-missing consensus base: `1`

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Samples: 278
- Consensus length: 124538
- Raw variant records considered: 2531
- Filtered SNP records available: 2022
- Filtered SNP sites applied inside track: 2022
- Raw-only failed variant sites masked: 508
- Missing consensus bases: 347258
- FASTA: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/11_callable_consensus/cpDNA.primary.callable_consensus.fa`
- Site table: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/11_callable_consensus/cpDNA.primary.callable_sites.tsv`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Samples: 278
- Consensus length: 44930
- Raw variant records considered: 192
- Filtered SNP records available: 146
- Filtered SNP sites applied inside track: 146
- Raw-only failed variant sites masked: 46
- Missing consensus bases: 152444
- FASTA: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/11_callable_consensus/mtDNA.primary.callable_consensus.fa`
- Site table: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/11_callable_consensus/mtDNA.primary.callable_sites.tsv`
