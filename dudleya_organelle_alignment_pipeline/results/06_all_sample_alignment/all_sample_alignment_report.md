# All-Sample Organelle Alignment

This step maps every primary paired-end sample to the combined
cpDNA/mtDNA reference and summarizes coverage using the Step 4
analysis tracks. It does not call variants or build consensus FASTAs.

## Inputs

- Sample table: `dudleya_organelle_alignment_pipeline/results/00_manifest/analysis_samples.tsv`
- Reference: `dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa`
- Analysis tracks: `dudleya_organelle_alignment_pipeline/results/05_analysis_masks/analysis_tracks.tsv`
- Minimum mapping quality retained in BAM/depth: `0`
- Minimum base quality used for depth: `20`

## Progress

- Samples completed in this run: 278
- Target samples for this invocation: 278

## Summary

- Total cpDNA+mtDNA mapped read records: 747507494
- Samples with initial QC notes: 3

## Median Breadth At 1x By Organelle

- Chloroplast: 1.000000
- Mitochondria: 0.998870

## Median Breadth At 1x By Analysis Track

- cpdna_duplicate_ir_mask: 1.000000
- cpdna_full_coverage: 1.000000
- cpdna_ir_regions: 1.000000
- cpdna_population_sites: 1.000000
- mtdna_high_confidence_unique: 1.000000
- mtdna_permissive_coverage: 0.998870

## Outputs

- `all_sample_alignment_sample_summary.tsv`: one row per sample.
- `all_sample_alignment_by_organelle.tsv`: one row per sample and organelle.
- `all_sample_alignment_by_track.tsv`: one row per sample and Step 4 track.
- `commands.tsv`: external commands run plus reuse decisions.
- `bam/`, `qc/`, and `logs/`: generated alignment artifacts ignored by git.

Review this report before variant calling, consensus generation, PCA,
tree building, Fst, or structure/admixture-style analyses.

## Samples With QC Notes

- `ABAB_MAD_LP_222_Du-589`: low_mitochondria_mapped_reads;low_chloroplast_breadth_ge_1x;low_mitochondria_breadth_ge_1x
- `CY_HUN_LP_265_Du-684`: low_mitochondria_breadth_ge_1x
- `CY_RED_LP_202_Du-561`: low_mitochondria_mapped_reads;low_chloroplast_breadth_ge_1x;low_mitochondria_breadth_ge_1x
