# Dudleya Organelle Alignment Preflight Summary

This is step 1 of the cpDNA/mtDNA alignment pipeline. It validates
sample naming, R1/R2 pairing, and population-code metadata before
any read alignment is attempted.

## Overall

- Samples discovered: 280
- Samples with exactly one R1 and one R2: 278
- Samples in primary paired-end alignment set: 278
- Samples excluded from primary paired-end alignment: 2
- Samples with resolved population metadata: 267
- Issues reported: 2

## Samples By Sequencing Batch

- QB3.Berkeley.241122/QB3.Dudleya.Results.250118: 8
- QB3.Berkeley.250811/QB3.250916.Results5genomes: 5
- QB3.Berkeley.251217/QB3.Results.260109: 267

## Samples By Naming Profile

- initial_du_dash: 8
- initial_du_lp: 5
- main_standard: 267

## R1/R2 Pair Status

- complete: 278
- missing_R1: 1
- missing_R2: 1

## Metadata Status

- resolved: 267
- unresolved_initial_sample: 13

## Primary Analysis Status

- exclude_missing_mate: 2
- include_primary_paired_end: 278

## Missing-Mate Policy

Samples without both mates are excluded from the primary paired-end
cpDNA/mtDNA alignment. They remain documented in `samples.tsv`,
`excluded_samples.tsv`, and `pairing_report.tsv`. If any missing-mate
sample is ever aligned as an individual single-end case, that run must
be reported separately as a sensitivity check and must not be mixed into
the primary paired-end dataset.

## Notes For The Next Pipeline Step

- `analysis_samples.tsv` is the input table for primary paired-end
  pilot read-to-reference alignment.
- `excluded_samples.tsv` records samples excluded from the primary
  alignment set and why.
- `main_standard` samples can be used for population-level analyses when
  their popcode appears in the population-code CSV.
- `initial_du_dash` and `initial_du_lp` samples can be aligned, but should
  remain metadata-unresolved until a manual lookup table is added.
- No alignment, trimming, variant calling, or consensus generation happens
  in this step.
