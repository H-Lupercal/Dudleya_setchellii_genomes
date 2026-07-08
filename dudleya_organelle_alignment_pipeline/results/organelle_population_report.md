# Dudleya Organelle Population Analysis Report

## Professor's Goal

"Create an alignment of all samples for the cpDNA and separately for the mtDNA
genomes using those annotated references; create a PCA visualization and a
phylogenetic tree (Neighbor Joining or maximum likelihood); generate population
genetic parameters such as admixture/structure plots with optimal K determined
empirically and Fst; think about whether to use Evan's pipeline for QC."

## Status

The requested primary organelle analysis has been generated for the 275-sample
downstream set. The workflow used the annotated cpDNA and mtDNA references from
`dudleya_organelle_reference_verification/`, excluded the two missing-mate
samples and three low-input QC-fail samples, and kept cpDNA and mtDNA outputs
separate.

## QC And Sample Set

- Biological sample rows discovered: 280.
- Complete paired-end samples initially alignable: 278.
- Missing-mate exclusions: 2.
- Downstream QC exclusions after all-sample mapping: 3.
- Primary downstream samples used for alignments, PCA, trees, and admixture: 275.
- Metadata-resolved populations used for Fst/population summaries: 34.

Evan Hackstadt's Snakemake pipeline was treated as the QC reference pattern:
sample-table driven processing, read-QC awareness, organized summaries, and
aggregate reporting. The core biological workflow here remains custom because
the analysis is organelle-specific, haploid, and uses cpDNA/mtDNA tracks rather
than a nuclear reference workflow.

## Alignments

Full callable-site consensus alignments:

| Organelle | Track | Samples | Sites | Output |
|---|---|---:|---:|---|
| cpDNA | `cpdna_population_sites` | 275 | 124,538 | `dudleya_organelle_alignment_pipeline/results/11_callable_consensus/cpDNA.primary.callable_consensus.fa` |
| mtDNA | `mtdna_high_confidence_unique` | 275 | 44,930 | `dudleya_organelle_alignment_pipeline/results/11_callable_consensus/mtDNA.primary.callable_consensus.fa` |

Filtered SNP-only alignments for PCA/admixture/Fst:

| Organelle | SNP Sites | Output |
|---|---:|---|
| cpDNA | 2,015 | `dudleya_organelle_alignment_pipeline/results/10_snp_alignment/cpDNA.primary.snp_alignment.fa` |
| mtDNA | 146 | `dudleya_organelle_alignment_pipeline/results/10_snp_alignment/mtDNA.primary.snp_alignment.fa` |

## PCA

| Organelle | Samples | SNPs | PC1 | PC2 | Plot |
|---|---:|---:|---:|---:|---|
| cpDNA | 275 | 2,015 | 36.62% | 14.65% | `dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.png` |
| mtDNA | 275 | 146 | 34.48% | 14.06% | `dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.png` |

Detailed PCA report:
`dudleya_organelle_alignment_pipeline/results/15_pca/primary.pca_report.md`

## Phylogenetic Trees

Maximum-likelihood trees were inferred with IQ-TREE using `GTR+F+G4` on the
full callable-site consensus alignments. The final tree run includes 1,000
ultrafast bootstrap replicates with `--bnni`.

| Organelle | Samples | Sites | Tree | Figure |
|---|---:|---:|---|---|
| cpDNA | 275 | 124,538 | `dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile` | `dudleya_organelle_alignment_pipeline/results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png` |
| mtDNA | 275 | 44,930 | `dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile` | `dudleya_organelle_alignment_pipeline/results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png` |

## Admixture-Style Clustering

ADMIXTURE was run separately for cpDNA and mtDNA across K=1..8 with
cross-validation and five seeded replicates per K. Because ADMIXTURE is
diploid-oriented, haploid organelle SNP calls were encoded as pseudo-diploid
homozygotes. Interpret these as organelle haplotype clustering plots, not
nuclear admixture.

| Organelle | Best K | Mean CV Error | CV SD | Structure Plot | CV Plot |
|---|---:|---:|---:|---|---|
| cpDNA | 8 | 0.08898600 | 0.01449154 | `dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png` | `dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/cpDNA.primary.admixture_cv.png` |
| mtDNA | 8 | 0.12644400 | 0.02207443 | `dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png` | `dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/mtDNA.primary.admixture_cv.png` |

Detailed admixture report:
`dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/primary.admixture_report.md`

## Fst And Population Summaries

Pairwise Fst and population summaries were generated for the 34 populations with
resolved population codes. The current Fst implementation uses a haploid
Nei-style differentiation estimate averaged across informative SNP sites.

| Organelle | Populations | Pairwise Comparisons | Pairwise Fst | Population Summary |
|---|---:|---:|---|---|
| cpDNA | 34 | 561 | `dudleya_organelle_alignment_pipeline/results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv` | `dudleya_organelle_alignment_pipeline/results/17_population_genetics/cpDNA.primary.population_genetics.population_summary.tsv` |
| mtDNA | 34 | 561 | `dudleya_organelle_alignment_pipeline/results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv` | `dudleya_organelle_alignment_pipeline/results/17_population_genetics/mtDNA.primary.population_genetics.population_summary.tsv` |

Detailed population-genetics report:
`dudleya_organelle_alignment_pipeline/results/17_population_genetics/primary.population_genetics_report.md`

## Tooling

Required analysis and visualization tools were installed or verified in the
local environment, including `bwa`, `samtools`, `bcftools`, `iqtree`, `plink`,
`admixture`, `vcftools`, `bedtools`, Python plotting/statistics packages, and R
plotting/tree packages. The current audit is:

`dudleya_organelle_alignment_pipeline/results/13_tool_audit/primary.tool_audit_report.md`

## Remaining Caveats

- Population summaries use only samples with resolved population codes, so
  initial unresolved DU-only samples are retained in alignments/PCA/trees but
  omitted from population-level Fst summaries.
- mtDNA results are limited to the high-confidence unique mtDNA track because
  whole-mtDNA mapping contains lower-confidence/repetitive regions.

## Method Notes

- ADMIXTURE is used here for organelle haplotype clustering, not nuclear
  admixture. Haploid organelle SNP calls are encoded as pseudo-diploid
  homozygotes for tool compatibility.
- The original fast ML topology trees remain in
  `dudleya_organelle_alignment_pipeline/results/12_phylogenetic_tree/`; the
  final tree deliverables are the bootstrap-supported trees in
  `dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/`.
