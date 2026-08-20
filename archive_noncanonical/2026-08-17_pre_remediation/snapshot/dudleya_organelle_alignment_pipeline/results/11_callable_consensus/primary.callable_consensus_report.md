# Callable-Site Consensus Alignment

This step builds full callable-site FASTA alignments for cpDNA and mtDNA.
Each alignment follows the Step 4 population-genetic BED track, starts
from the annotated organelle reference, overlays Step 8 filtered haploid
SNP genotypes, masks Step 7 raw variant sites that failed filtering, and
uses Step 5 depth files to write `N` at bases below the minimum depth.

## Run

- Run label: `primary`
- Minimum depth for a non-missing consensus base: `1`

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Samples: 275
- Consensus length: 124538
- Raw variant records considered: 2475
- Filtered SNP records available: 2015
- Filtered SNP sites applied inside track: 2015
- Raw-only failed variant sites masked: 459
- Missing consensus bases: 127485
- FASTA: `dudleya_organelle_alignment_pipeline/results/11_callable_consensus/cpDNA.primary.callable_consensus.fa`
- Site table: `dudleya_organelle_alignment_pipeline/results/11_callable_consensus/cpDNA.primary.callable_sites.tsv`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Samples: 275
- Consensus length: 44930
- Raw variant records considered: 190
- Filtered SNP records available: 146
- Filtered SNP sites applied inside track: 146
- Raw-only failed variant sites masked: 44
- Missing consensus bases: 31313
- FASTA: `dudleya_organelle_alignment_pipeline/results/11_callable_consensus/mtDNA.primary.callable_consensus.fa`
- Site table: `dudleya_organelle_alignment_pipeline/results/11_callable_consensus/mtDNA.primary.callable_sites.tsv`
