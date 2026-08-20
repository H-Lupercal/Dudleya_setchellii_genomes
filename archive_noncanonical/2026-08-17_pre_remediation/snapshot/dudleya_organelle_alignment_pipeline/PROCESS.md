# Analysis Process And Stage Index

This is the authoritative, ordered map of the cpDNA/mtDNA analysis and the single
source of truth for stage order, inputs, and outputs.

Each stage is identified by its `results/NN_.../` directory number, which is the
canonical stage identifier used in all paths and documentation. Directory numbers
are contiguous (00-21) and include the pilot investigations (03, 04) and the
separate initial/final runs. The auto-generated stage reports use descriptive
titles (for example `# Phylogenetic Trees`).

The completed analysis was run under the `primary` run label on the 275-sample
downstream set. A 5-sample `smoke` run label exists for controlled validation of
the variant-calling stage and is not part of the primary results.

## Stage Map

| Dir | Purpose | Runner script | Output directory | Report / summary |
|---|---|---|---|---|
| 00 | Sample manifest and preflight validation | `scripts/build_sample_manifest.py` | `results/00_manifest/` | `preflight_summary.md` |
| 01 | Reference validation, indexing, and pilot-sample selection | `scripts/prepare_reference_and_pilot.py` | `results/01_reference_pilot/` | `reference_pilot_summary.md` |
| 02 | Pilot organelle alignment (15 samples) | `scripts/run_pilot_alignment.py` | `results/02_pilot_alignment/` | `pilot_alignment_report.md` |
| 03 | mtDNA repeat/placement investigation (pilot phase) | investigation evidence; consumed by Stage 05 | `results/03_mtdna_investigation/` | `mtdna_investigation_report.md` |
| 04 | cpDNA inverted-repeat verification (pilot phase) | investigation evidence; consumed by Stage 05 | `results/04_cpdna_investigation/` | `cpdna_verification_report.md` |
| 05 | Analysis masks and coverage/population tracks | `scripts/build_analysis_masks.py` | `results/05_analysis_masks/` | `mask_summary.md` |
| 06 | All-sample organelle alignment and track-aware QC | `scripts/run_all_sample_alignment.py` | `results/06_all_sample_alignment/` | `all_sample_alignment_report.md` |
| 07 | Downstream include/exclude sample set (275 samples) | `scripts/build_downstream_sample_set.py` | `results/07_downstream_sample_set/` | `downstream_sample_set_report.md` |
| 08 | Haploid variant calling (cpDNA and mtDNA separately) | `scripts/run_variant_calling.py` | `results/08_variant_calling/` | `primary.variant_calling_report.md` |
| 09 | Variant filtering and callable-site definition | `scripts/run_variant_filtering.py` | `results/09_variant_filtering/` | `primary.variant_filtering_report.md` |
| 10 | Filtered haploid SNP alignments | `scripts/run_snp_alignment.py` | `results/10_snp_alignment/` | `primary.snp_alignment_report.md` |
| 11 | Full callable-site consensus alignments | `scripts/run_callable_consensus.py` | `results/11_callable_consensus/` | `primary.callable_consensus_report.md` |
| 12 | Maximum-likelihood trees, initial (no bootstrap) | `scripts/run_phylogenetic_tree.py` | `results/12_phylogenetic_tree/` | `primary.phylogenetic_tree_report.md` |
| 13 | Bioinformatics tool audit | `scripts/run_tool_audit.py` | `results/13_tool_audit/` | `primary.tool_audit_report.md` |
| 14 | Tree figures for the Stage 12 trees | `scripts/run_tree_visualization.py` | `results/14_tree_visualization/` | `primary.tree_visualization_report.md` |
| 15 | PCA (cpDNA and mtDNA) | `scripts/run_pca_analysis.py` | `results/15_pca/` | `primary.pca_report.md` |
| 16 | ADMIXTURE-style clustering, single-run | `scripts/run_admixture_analysis.py` | `results/16_admixture/` | `primary.admixture_report.md` |
| 17 | Pairwise Fst and population summaries | `scripts/run_population_genetics.py` | `results/17_population_genetics/` | `primary.population_genetics_report.md` |
| 18 | ADMIXTURE-style clustering, five-replicate (final) | `scripts/run_admixture_analysis.py --replicates 5 --output-dir results/18_admixture_replicates` | `results/18_admixture_replicates/` | `primary.admixture_report.md` |
| 19 | Maximum-likelihood trees, 1,000 UFBoot (final) | `scripts/run_phylogenetic_tree.py --bootstrap-replicates 1000 --output-dir results/19_bootstrap_phylogenetic_tree` | `results/19_bootstrap_phylogenetic_tree/` | `primary.phylogenetic_tree_report.md` |
| 20 | Tree figures for the Stage 19 bootstrap trees (final) | `scripts/run_tree_visualization.py` (pointed at the Stage 19 trees) | `results/20_bootstrap_tree_visualization/` | `primary.tree_visualization_report.md` |
| 21 | Haploid cpDNA/mtDNA haplotype networks and PopART exports | `scripts/run_haplotype_network.py` | `results/21_haplotype_network/` | `primary.haplotype_network_report.md` |

## Initial-Versus-Final Runs

Three analyses were run twice: a first pass and a final, more rigorous pass. The
final pass supersedes the first for interpretation and for the deliverables.

| First pass | Final pass | Difference |
|---|---|---|
| Stage 12 (ML trees) | Stage 19 (ML trees) | Final adds 1,000 ultrafast bootstrap replicates with BNNI. |
| Stage 14 (tree figures) | Stage 20 (tree figures) | Final renders the bootstrap-supported Stage 19 trees. |
| Stage 16 (ADMIXTURE) | Stage 18 (ADMIXTURE) | Final uses five seeded replicates per K and selects K by lowest mean CV error. |

## Final Deliverables

The final deliverables and their exact paths are listed in
`results/final_deliverables_manifest.tsv`, and the integrated methods/results
narrative is `results/organelle_population_report.md`.

Stage 21 is additive: it provides a haploid-native visualization alongside the
existing Stage 18 ADMIXTURE-style results rather than replacing them.

## Command Provenance

Each stage has a runner script and a `report.md` (or `*_summary.md`) recording
its parameters and results. Stages that shell out to external tools also write a
verbatim `commands.tsv`: mapping (02, 06), variant calling (08), variant
filtering (09), trees (12, 19), and haplotype networks (21). For the remaining stages, the exact
invocation is given in the reproduction command blocks in `README.md` and the
top-level `../README.md`, and the key parameters are recorded in each stage
report.
