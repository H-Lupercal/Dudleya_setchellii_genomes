#!/usr/bin/env bash
set -euo pipefail

THREADS="${1:-16}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_DIR="$ROOT/full_pipeline_run"
RESULTS="$RUN_DIR/results"
LOG_DIR="$RUN_DIR/logs"
TOOL_PATH="$ROOT/.tools/bioconda-env/bin:$PATH"
PYTHON="${PYTHON:-$ROOT/.tools/bioconda-env/bin/python3}"

mkdir -p "$RESULTS" "$LOG_DIR"

run_stage() {
  local stage_id="$1"
  local description="$2"
  shift 2
  local log_file="$LOG_DIR/${stage_id}.log"
  {
    printf '=== %s %s ===\n' "$stage_id" "$description"
    printf 'Started: %s\n' "$(date -Is)"
    printf 'Command:'
    printf ' %q' "$@"
    printf '\n\n'
    "$@"
    printf '\nFinished: %s\n' "$(date -Is)"
  } >"$log_file" 2>&1
  printf '%s\t%s\t%s\n' "$stage_id" "$description" "$log_file" >>"$LOG_DIR/stage_status.tsv"
}

cd "$ROOT"
{
  printf '\nResumed from Stage 07: %s\n' "$(date -Is)"
  printf 'Threads: %s\n' "$THREADS"
} >>"$RUN_DIR/run_metadata.txt"

run_stage 07_downstream_sample_set "Downstream include/exclude sample set" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/build_downstream_sample_set.py \
  --analysis-sample-table "$RESULTS/00_manifest/analysis_samples.tsv" \
  --upstream-excluded-table "$RESULTS/00_manifest/excluded_samples.tsv" \
  --qc-decision-table "$RESULTS/06_all_sample_alignment/all_sample_alignment_sample_summary.tsv" \
  --output-dir "$RESULTS/07_downstream_sample_set" \
  --expected-included -1

PATH="$TOOL_PATH" run_stage 08_variant_calling "Haploid cpDNA and mtDNA variant calling" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_variant_calling.py \
  --sample-table "$RESULTS/07_downstream_sample_set/included_samples.tsv" \
  --track-table "$RESULTS/05_analysis_masks/analysis_tracks.tsv" \
  --bam-dir "$RESULTS/06_all_sample_alignment/bam" \
  --output-dir "$RESULTS/08_variant_calling" \
  --run-label primary \
  --threads "$THREADS" \
  --force

PATH="$TOOL_PATH" run_stage 09_variant_filtering "Variant filtering and callable-site definition" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_variant_filtering.py \
  --variant-calling-dir "$RESULTS/08_variant_calling" \
  --output-dir "$RESULTS/09_variant_filtering" \
  --run-label primary \
  --threads "$THREADS" \
  --force

run_stage 10_snp_alignment "Filtered haploid SNP alignments" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_snp_alignment.py \
  --variant-filtering-dir "$RESULTS/09_variant_filtering" \
  --output-dir "$RESULTS/10_snp_alignment" \
  --run-label primary

run_stage 11_callable_consensus "Full callable-site consensus alignments" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_callable_consensus.py \
  --sample-table "$RESULTS/07_downstream_sample_set/included_samples.tsv" \
  --track-table "$RESULTS/05_analysis_masks/analysis_tracks.tsv" \
  --depth-dir "$RESULTS/06_all_sample_alignment/qc" \
  --variant-calling-dir "$RESULTS/08_variant_calling" \
  --variant-filtering-dir "$RESULTS/09_variant_filtering" \
  --output-dir "$RESULTS/11_callable_consensus" \
  --run-label primary

PATH="$TOOL_PATH" run_stage 12_phylogenetic_tree "Initial maximum-likelihood trees" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_phylogenetic_tree.py \
  --consensus-dir "$RESULTS/11_callable_consensus" \
  --output-dir "$RESULTS/12_phylogenetic_tree" \
  --run-label primary \
  --threads "$THREADS" \
  --force

PATH="$TOOL_PATH" run_stage 13_tool_audit "Bioinformatics tool audit" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_tool_audit.py \
  --output-dir "$RESULTS/13_tool_audit" \
  --audit-label primary

run_stage 14_tree_visualization "Tree figures for initial ML trees" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_tree_visualization.py \
  --tree-dir "$RESULTS/12_phylogenetic_tree" \
  --output-dir "$RESULTS/14_tree_visualization" \
  --run-label primary

run_stage 15_pca "PCA visualizations" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_pca_analysis.py \
  --snp-alignment-dir "$RESULTS/10_snp_alignment" \
  --metadata-path "$RESULTS/07_downstream_sample_set/included_samples.tsv" \
  --output-dir "$RESULTS/15_pca" \
  --run-label primary

PATH="$TOOL_PATH" run_stage 16_admixture "Single-run ADMIXTURE-style clustering" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_admixture_analysis.py \
  --snp-alignment-dir "$RESULTS/10_snp_alignment" \
  --metadata-path "$RESULTS/07_downstream_sample_set/included_samples.tsv" \
  --output-dir "$RESULTS/16_admixture" \
  --run-label primary \
  --threads "$THREADS" \
  --force

run_stage 17_population_genetics "Pairwise Fst and population summaries" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_population_genetics.py \
  --snp-alignment-dir "$RESULTS/10_snp_alignment" \
  --metadata-path "$RESULTS/07_downstream_sample_set/included_samples.tsv" \
  --output-dir "$RESULTS/17_population_genetics" \
  --run-label primary

PATH="$TOOL_PATH" run_stage 18_admixture_replicates "Five-replicate final ADMIXTURE-style clustering" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_admixture_analysis.py \
  --snp-alignment-dir "$RESULTS/10_snp_alignment" \
  --metadata-path "$RESULTS/07_downstream_sample_set/included_samples.tsv" \
  --output-dir "$RESULTS/18_admixture_replicates" \
  --run-label primary \
  --threads "$THREADS" \
  --replicates 5 \
  --force

PATH="$TOOL_PATH" run_stage 19_bootstrap_phylogenetic_tree "Final 1,000-UFBoot maximum-likelihood trees" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_phylogenetic_tree.py \
  --consensus-dir "$RESULTS/11_callable_consensus" \
  --output-dir "$RESULTS/19_bootstrap_phylogenetic_tree" \
  --run-label primary \
  --threads "$THREADS" \
  --bootstrap-replicates 1000 \
  --force

run_stage 20_bootstrap_tree_visualization "Tree figures for final bootstrap ML trees" \
  "$PYTHON" dudleya_organelle_alignment_pipeline/scripts/run_tree_visualization.py \
  --tree-dir "$RESULTS/19_bootstrap_phylogenetic_tree" \
  --output-dir "$RESULTS/20_bootstrap_tree_visualization" \
  --run-label primary

{
  printf '\nFinished: %s\n' "$(date -Is)"
  printf 'Result size: '
  du -sh "$RUN_DIR" | awk '{print $1}'
} >>"$RUN_DIR/run_metadata.txt"
