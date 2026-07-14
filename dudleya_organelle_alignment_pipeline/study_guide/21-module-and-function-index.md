# Chapter 21 — Module and Function Index

> Part 4 of 4 · Practice and Reference · Prev:
> [Glossary](./20-glossary.md) · Next: [External Tool and Format
> Reference](./22-external-tool-and-format-reference.md)

This is the coverage audit: every implementation module, every runner script, and
every top-level function, so you can confirm nothing in the pipeline is
undocumented and locate any symbol fast. **Scaffolding** functions repeat across
modules and are taught once in [Chapter 3](./03-reusable-code-patterns.md):
`read_tsv`, `write_tsv` / `write_dataclass_tsv`, `open_text`, `build_arg_parser`,
`main`. They are listed per module for completeness but not re-explained.

## Runner scripts (`scripts/`)

All sixteen are identical thin wrappers ([Chapter 3, §3.1](./03-reusable-code-patterns.md))
that fix `sys.path` and call their module's `main`: `build_sample_manifest.py`,
`prepare_reference_and_pilot.py`, `run_pilot_alignment.py`,
`build_analysis_masks.py`, `run_all_sample_alignment.py`,
`build_downstream_sample_set.py`, `run_variant_calling.py`,
`run_variant_filtering.py`, `run_snp_alignment.py`, `run_callable_consensus.py`,
`run_phylogenetic_tree.py`, `run_tool_audit.py`, `run_tree_visualization.py`,
`run_pca_analysis.py`, `run_admixture_analysis.py`, `run_population_genetics.py`.

## `__init__.py`

Package marker; docstring only, no code.

## `manifest.py` — Stage 00 · [Chapter 7](./07-manifest-and-reference-preflight.md)

Classes: `PopulationCode`, `FastqRecord`, `ManifestRow`, `ManifestIssue`.
Functions: `clean_cell`, `discover_fastq_paths`, `infer_batch`, `classify_prefix`,
`parse_fastq_path`, `load_population_codes`, `infer_species_from_popcode`,
`build_manifest`, `determine_pair_status`, `determine_analysis_status`,
`join_unique`, `join_paths`, `write_manifest_outputs`, `primary_analysis_rows`,
`excluded_analysis_rows`, `write_dataclass_tsv`, `write_preflight_summary`,
`format_counts`, `build_arg_parser`, `main`.

## `prepare_reference_and_pilot.py` — Stage 01 · [Chapter 7](./07-manifest-and-reference-preflight.md)

Classes: `ReferenceValidationError`, `ReferenceCheck`, `ToolCheck`, `IndexCheck`.
Functions: `read_fasta_lengths`, `validate_reference_records`, `check_tools`,
`prepare_indexes`, `run_command`, `bwa_suffixes`, `read_tsv`,
`select_pilot_samples`, `main_standard_candidates_by_species`, `with_reason`,
`species_group`, `sample_sort_key`, `population_sort_key`, `safe_reason`,
`write_pilot_samples`, `write_dataclass_tsv`, `write_summary`, `build_arg_parser`,
`main`.

## `pilot_alignment.py` — Stage 02 · [Chapter 8](./08-pilot-mapping-and-investigations.md)

The definitional home of the alignment/QC helpers reused by Stages 06, 08, 11.
Classes: `AlignmentError`, `AlignmentSample`, `AlignmentOutputs`,
`OrganelleMetrics`. Functions: `safe_sample_name`, `read_tsv`, `split_path_field`,
`read_alignment_samples`, `count_fastq_records`, `read_fai_lengths`,
`parse_idxstats_file`, `parse_depth_file`, `fmt_float`, `build_sample_summary`,
`build_organelle_summary_rows`, `outputs_for_sample`, `outputs_are_ready`,
`require_tools`, `require_reference_indexes`, `build_depth_command`, `shlex_join`,
`run_alignment_commands`, `run_qc_commands`, `run_logged_command`, `write_tsv`,
`write_report`, `median_breadth`, `run_pilot_alignment`, `build_arg_parser`,
`main`.

## `analysis_masks.py` — Stage 05 · [Chapter 9](./09-masks-alignment-and-sample-qc.md)

Classes: `MaskDefinitionError`, `Region`, `CpdnaTracks`, `MtdnaTracks`.
Functions: `interval_to_bed_fields`, `validate_interval`,
`validate_interval_within_reference`, `read_tsv`, `build_cpdna_tracks`,
`read_major_cpdna_repeat_pair`, `build_mtdna_tracks`,
`read_mtdna_high_confidence_regions`, `complement_regions`, `build_named_region`,
`merge_coordinate_pairs`, `generate_analysis_masks`, `write_bed`,
`write_region_manifest`, `write_track_manifest`, `write_summary`,
`build_arg_parser`, `main`.

## `all_sample_alignment.py` — Stage 06 · [Chapter 9](./09-masks-alignment-and-sample-qc.md)

Imports the alignment core from `pilot_alignment.py`. Classes: `TrackRegion`,
`TrackMetrics`. Functions: `read_tsv`, `read_track_regions`,
`validate_track_regions`, `parse_track_depth_file`, `initialize_track_counters`,
`build_track_summary_rows`, `track_order`, `run_all_sample_alignment`,
`write_all_sample_outputs`, `write_report`, `median_breadth`,
`format_track_medians`, `build_arg_parser`, `main`.

## `downstream_sample_set.py` — Stage 07 · [Chapter 9](./09-masks-alignment-and-sample-qc.md)

Classes: `DownstreamSampleSetError`. Functions: `read_tsv`, `write_tsv`,
`build_downstream_sample_set`, `build_included_row`, `build_qc_excluded_row`,
`build_upstream_excluded_row`, `append_evidence`,
`write_downstream_sample_set_outputs`, `write_report`,
`generate_downstream_sample_set`, `build_arg_parser`, `main`.

## `variant_calling.py` — Stage 08 · [Chapter 10](./10-variants-to-alignments.md)

Home of `labeled_output_name`, reused by every downstream stage. Classes:
`VariantCallingError`, `VariantSample`, `VariantTrack`, `VariantCallInputs`,
`VariantCallResult`. Functions: `read_tsv`, `write_tsv`, `read_variant_samples`,
`read_variant_tracks`, `label_output_prefix`, `write_variant_call_inputs`,
`build_bcftools_commands`, `require_bcftools`, `run_variant_call_for_track`,
`outputs_are_ready`, `count_vcf_records`, `run_variant_calling`,
`labeled_output_name`, `write_variant_calling_outputs`, `write_report`,
`build_arg_parser`, `main`.

## `variant_filtering.py` — Stage 09 · [Chapter 10](./10-variants-to-alignments.md)

Classes: `VariantFilteringError`, `FilterInput`, `FilterResult`. Functions:
`read_tsv`, `write_tsv`, `require_bcftools`, `read_filter_inputs`,
`output_filtered_vcf_path`, `output_log_path`, `build_filter_command`,
`build_index_command`, `outputs_are_ready`, `filter_one_input`,
`run_variant_filtering`, `write_filtering_outputs`, `write_report`,
`build_arg_parser`, `main`.

## `snp_alignment.py` — Stage 10 · [Chapter 10](./10-variants-to-alignments.md)

Classes: `SnpAlignmentError`, `SnpAlignment`, `SnpAlignmentInput`,
`SnpAlignmentResult`. Functions: `read_tsv`, `write_tsv`, `open_text`,
`read_alignment_inputs`, `genotype_to_base`, `build_snp_alignment`, `write_fasta`,
`write_site_table`, `alignment_output_paths`, `build_one_alignment`,
`run_snp_alignment`, `write_alignment_outputs`, `write_report`, `build_arg_parser`,
`main`.

## `callable_consensus.py` — Stage 11 · [Chapter 10](./10-variants-to-alignments.md)

Classes: `CallableConsensusError`, `BedInterval`, `VcfVariant`,
`CallableAlignment`, `ConsensusInput`, `ConsensusResult`. Functions: `read_tsv`,
`write_tsv`, `open_text`, `read_fasta`, `read_bed`, `read_sample_names`,
`read_vcf_variants`, `genotype_to_base`, `build_track_template`,
`read_covered_indexes`, `depth_path_for_sample`, `build_callable_consensus`,
`read_consensus_inputs`, `write_fasta`, `write_site_table`,
`consensus_output_paths`, `build_one_consensus`, `run_callable_consensus`,
`write_consensus_outputs`, `write_report`, `build_arg_parser`, `main`.

## `phylogenetic_tree.py` — Stages 12/19 · [Chapter 11](./11-phylogenetic-trees.md)

Classes: `PhylogeneticTreeError`, `TreeInput`, `TreeResult`. Functions: `read_tsv`,
`write_tsv`, `require_iqtree`, `read_tree_inputs`, `tree_output_prefix`,
`build_iqtree_command`, `tree_outputs_ready`, `run_one_tree`,
`run_phylogenetic_trees`, `write_tree_outputs`, `write_report`, `build_arg_parser`,
`main`.

## `tree_visualization.py` — Stages 14/20 · [Chapter 11](./11-phylogenetic-trees.md)

Classes: `TreeVisualizationError`, `TreeFigureInput`, `TreeFigureResult`.
Functions: `read_tsv`, `write_tsv`, `read_tree_figure_inputs`,
`compute_tree_figure_size`, `render_tree_figure`,
`write_tree_visualization_outputs`, `write_report`, `run_tree_visualizations`,
`build_arg_parser`, `main`.

## `pca_analysis.py` — Stage 15 · [Chapter 12](./12-pca-and-clustering.md)

Home of `read_fasta`, `read_sample_metadata`, `choose_plot_group`, reused by
Stages 16/18 and 17. Classes: `PcaAnalysisError`, `PcaInput`, `PcaResult`.
Functions: `read_tsv`, `write_tsv`, `read_pca_inputs`, `read_fasta`,
`build_haploid_snp_matrix`, `read_sample_metadata`, `choose_plot_group`,
`run_pca`, `write_pca_tables`, `write_pca_plot`, `run_one_pca`, `write_pca_outputs`,
`write_pca_report`, `run_pca_analysis`, `build_arg_parser`, `main`.

## `admixture_analysis.py` — Stages 16/18 · [Chapter 12](./12-pca-and-clustering.md)

Classes: `AdmixtureAnalysisError`, `AdmixtureInput`. Functions: `read_tsv`,
`write_tsv`, `require_admixture`, `require_plink`, `read_admixture_inputs`,
`write_pseudo_diploid_ped_map`, `build_admixture_command`,
`build_plink_make_bed_command`, `run_plink_make_bed`, `parse_cv_error`,
`summarize_replicate_stability`, `run_admixture_for_k`, `read_q_matrix`,
`write_q_table_and_plot`, `write_cv_plot`, `write_admixture_outputs`,
`write_admixture_report`, `run_admixture_analysis`, `build_arg_parser`, `main`.

## `population_genetics.py` — Stage 17 · [Chapter 13](./13-population-fst.md)

Classes: `PopulationGeneticsError`, `PopulationInput`, `PopulationResult`.
Functions: `read_tsv`, `write_tsv`, `read_population_inputs`,
`population_code_for_sample`, `group_sequences_by_population`,
`allele_counts_at_site`, `gene_diversity`, `compute_pairwise_fst`,
`compute_haplotype_diversity`, `compute_nucleotide_diversity`,
`private_variant_count`, `run_one_population_summary`, `run_population_genetics`,
`write_population_genetics_outputs`, `write_population_genetics_report`,
`build_arg_parser`, `main`.

## `tool_audit.py` — Stage 13 · [Chapter 14](./14-tool-audit.md)

Classes: `ToolSpec`, `ToolResult`, `AuditSummary`. Module constant: `TOOL_SPECS`
(the 24-tool audit list). Functions: `run_version_command`, `first_version_line`,
`missing_note`, `check_tool`, `run_tool_audit`, `summarize_audit`, `write_tsv`,
`write_tool_audit_outputs`, `write_report`, `build_arg_parser`, `main`.

## Test files (`tests/`)

One per module, each a worked example ([Chapter 3, §3.8](./03-reusable-code-patterns.md)):
`test_manifest.py`, `test_prepare_reference_and_pilot.py`, `test_pilot_alignment.py`,
`test_analysis_masks.py`, `test_all_sample_alignment.py`,
`test_downstream_sample_set.py`, `test_variant_calling.py`,
`test_variant_filtering.py`, `test_snp_alignment.py`, `test_callable_consensus.py`,
`test_phylogenetic_tree.py`, `test_tree_visualization.py`, `test_pca_analysis.py`,
`test_admixture_analysis.py`, `test_population_genetics.py`, `test_tool_audit.py`.

## Symbol → chapter quick lookup

| To find… | Look in |
|---|---|
| filename parsing (`classify_prefix`, regexes) | [Ch 7](./07-manifest-and-reference-preflight.md) |
| the `bwa/samtools` pipe (`run_alignment_commands`) | [Ch 8](./08-pilot-mapping-and-investigations.md) |
| coordinate conversion (`interval_to_bed_fields`) | [Ch 5](./05-bioinformatics-file-formats.md), [Ch 9](./09-masks-alignment-and-sample-qc.md) |
| IR masking (`build_cpdna_tracks`, `complement_regions`) | [Ch 9](./09-masks-alignment-and-sample-qc.md) |
| the `bcftools` calls (`build_bcftools_commands`, `build_filter_command`) | [Ch 10](./10-variants-to-alignments.md), [Ch 4](./04-shell-and-external-tools.md) |
| genotype decoding (`genotype_to_base`) | [Ch 5](./05-bioinformatics-file-formats.md), [Ch 10](./10-variants-to-alignments.md) |
| callable consensus (`build_callable_consensus`) | [Ch 10](./10-variants-to-alignments.md), [Ch 23](./23-capstone-sample-trace.md) |
| IQ-TREE command (`build_iqtree_command`) | [Ch 11](./11-phylogenetic-trees.md) |
| PCA matrix (`build_haploid_snp_matrix`) | [Ch 12](./12-pca-and-clustering.md) |
| pseudo-diploid encoding (`write_pseudo_diploid_ped_map`) | [Ch 12](./12-pca-and-clustering.md) |
| K selection (`summarize_replicate_stability`) | [Ch 12](./12-pca-and-clustering.md) |
| Fst (`compute_pairwise_fst`, `gene_diversity`) | [Ch 13](./13-population-fst.md) |
| dependency injection (`check_tool`) | [Ch 14](./14-tool-audit.md) |
| run-label naming (`labeled_output_name`) | [Ch 3](./03-reusable-code-patterns.md) |

> Next: [Chapter 22 — External Tool and File-Format Reference](./22-external-tool-and-format-reference.md)
