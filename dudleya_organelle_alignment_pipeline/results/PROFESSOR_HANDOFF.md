# Professor Handoff

Start here:

`dudleya_organelle_alignment_pipeline/results/organelle_population_report.md`

That report states the professor's requested goal, summarizes the sample/QC
decisions, and links the cpDNA and mtDNA alignments, PCA plots, ML trees,
admixture/structure-style plots, Fst tables, population summaries, and tool
audit.

The machine-readable file index is:

`dudleya_organelle_alignment_pipeline/results/final_deliverables_manifest.tsv`

Primary headline results:

- 275 downstream samples were used for cpDNA/mtDNA alignments, PCA, trees, and
  ADMIXTURE-style clustering.
- cpDNA callable alignment: 124,538 sites.
- mtDNA callable alignment: 44,930 sites.
- cpDNA PCA: 2,015 SNPs; PC1 36.62%; PC2 14.65%.
- mtDNA PCA: 146 SNPs; PC1 34.48%; PC2 14.06%.
- cpDNA ML tree and mtDNA ML tree were inferred with IQ-TREE and 1,000
  ultrafast bootstrap replicates.
- cpDNA ADMIXTURE-style best K: 8.
- mtDNA ADMIXTURE-style best K: 8.
- Fst/population summaries cover 34 metadata-resolved populations and 561
  pairwise population comparisons per organelle.

Important caveats to preserve when sharing:

- Fst/population summaries include only samples with resolved population codes.
- mtDNA population genetics uses the high-confidence unique mtDNA track.

Method notes:

- ADMIXTURE was run as organelle haplotype clustering, not nuclear admixture.
  Haploid calls were encoded as pseudo-diploid homozygotes.
- ADMIXTURE K selection now uses five seeded replicates per K and selects K by
  lowest mean cross-validation error.
