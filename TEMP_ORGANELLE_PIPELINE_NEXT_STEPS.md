# Temporary Organelle Pipeline Next Steps

This temporary queue is now intentionally empty. The primary organelle analysis
requested by the professor has been executed through the final handoff outputs.

Use these files as the current source of truth:

```text
dudleya_organelle_alignment_pipeline/results/PROFESSOR_HANDOFF.md
dudleya_organelle_alignment_pipeline/results/organelle_population_report.md
dudleya_organelle_alignment_pipeline/results/final_deliverables_manifest.tsv
```

Completed primary outputs:

- cpDNA and mtDNA callable consensus alignments.
- cpDNA and mtDNA SNP alignments.
- cpDNA and mtDNA PCA plots.
- cpDNA and mtDNA maximum-likelihood trees with 1,000 ultrafast bootstraps.
- Bootstrap tree figures.
- ADMIXTURE-style organelle haplotype clustering with five seeded replicates
  per K.
- Pairwise Fst and per-population summaries.
- Tool audit and final professor-facing reports.

Known interpretation notes:

- ML trees were delivered; no separate rendered NJ quick-check figure was made.
- Fst is a custom haploid Nei-style estimate over informative SNPs.
- Fst/population summaries include only metadata-resolved populations.
- mtDNA popgen uses the high-confidence unique mtDNA track, not the full
  repeat-rich mitochondrial contig.

Add new tasks below only if new follow-up work is requested.
