# Pilot Organelle Alignment

This step aligns the representative pilot samples to the combined
cpDNA/mtDNA reference and summarizes organelle mapping signal. It does
not call variants, make consensus FASTAs, or build final alignments.

## Inputs

- Reference: `dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa`
- Minimum mapping quality retained in BAM/depth: `0`
- Minimum base quality used for depth: `20`

## Summary

- Samples attempted: 15
- Total cpDNA+mtDNA mapped read records: 46584669
- Samples with initial QC notes: 1

## Median Breadth At 1x

- Chloroplast: 0.999993
- Mitochondria: 0.960445

## Outputs

- `pilot_alignment_sample_summary.tsv`: one row per sample.
- `pilot_alignment_by_organelle.tsv`: one row per sample and organelle.
- `commands.tsv`: external commands run plus any reuse decisions.
- `bam/`: filtered, sorted, indexed organelle BAM files.
- `qc/`: per-sample flagstat, idxstats, and depth files.

Review the sample and organelle summaries before scaling to all primary
paired-end samples.

## Samples With QC Notes

- `ABAB_MAD_LP_222_Du-589`: low_mitochondria_mapped_reads;low_chloroplast_breadth_ge_1x;low_mitochondria_breadth_ge_1x
