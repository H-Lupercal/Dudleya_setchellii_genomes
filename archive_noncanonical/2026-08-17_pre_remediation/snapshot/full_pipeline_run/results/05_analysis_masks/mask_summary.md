# Analysis Masks

This step defines the cpDNA and mtDNA tracks that the all-sample run
must use. It does not align reads, call variants, or create final
population-genetic outputs.

## Coordinate Systems

- BED files are 0-based, half-open.
- `analysis_regions.tsv` records the same intervals as 1-based inclusive
  coordinates plus their BED coordinates.

## cpDNA Tracks

- Source: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/04_cpdna_investigation/cpdna_self_repeat_intervals.tsv`
- Full cpDNA reference length: 150274 bp.
- Strategy: keep one chloroplast IR copy for population-genetic outputs.
- Duplicate IR bases masked from cpDNA population sites: 25736.
- cpDNA population-site bases retained: 124538.
- Use `cpdna_full_coverage_regions.bed` for sample-level coverage QC.
- Use `cpdna_population_sites.bed` for PCA, Fst, trees, and
  admixture-style clustering inputs.

## mtDNA Tracks

- Source: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/03_mtdna_investigation/mtdna_high_mapq_consensus_intervals.tsv`
- Full mtDNA reference length: 243359 bp.
- Strategy: keep mtDNA in two tracks.
- Use `mtdna_permissive_coverage_regions.bed` for sample-level
  permissive MAPQ coverage QC.
- Use `mtdna_high_confidence_unique_regions.bed` for mtDNA variant
  calling and population genetics.
- High-confidence mtDNA threshold: intervals supported by at least
  12 usable pilot samples.
- High-confidence mtDNA bases retained: 44930.

## Outputs

- `analysis_tracks.tsv`: machine-readable track purpose and downstream use.
- `analysis_regions.tsv`: machine-readable interval audit table.
- `*.bed`: regions and masks consumed by later alignment/QC/variant steps.
