#!/usr/bin/env python3
"""Build canonical reports, invalidation accounting, manifests, and acceptance gate."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
import tomllib
from dataclasses import asdict
from pathlib import Path

from Bio import Phylo
from organelle_pipeline.configuration import validate_publication_config
from organelle_pipeline.inventory import ACCEPTABLE_SOURCE_VALIDATION_STATUSES, validate_inventory
from organelle_pipeline.paths import repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    pipeline_code_digest,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)
from organelle_pipeline.references import read_single_fasta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def invalidation_for(path: str, run_id: str) -> tuple[str, str, str]:
    lower = path.lower()
    base = "canonical_publication"
    if lower.endswith((".py", ".sh", ".r", ".yml", ".yaml")):
        return (
            "superseded_code",
            "legacy implementation retained for audit but unsupported",
            f"{base}/pipeline/",
        )
    if "population_genetics" in lower:
        return (
            "invalidated",
            "SNP-only nucleotide-diversity denominator and nonstandard clamped FST",
            f"{base}/results/popgen/{run_id}/",
        )
    if "callable_consensus" in lower or "snp_alignment" in lower:
        return (
            "invalidated",
            "consensus inherited permissive mapping and reference-filled uncertain sites",
            f"{base}/results/alignments/{run_id}/",
        )
    if "variant" in lower:
        return (
            "invalidated",
            "legacy calls lacked canonical GQ/DP genotype masking and singleton policy",
            f"{base}/results/variants/{run_id}/",
        )
    if "admixture" in lower:
        return (
            "retired_to_supplement",
            "linked haploid organelle markers violate primary ADMIXTURE assumptions",
            f"{base}/results/supplement/{run_id}/admixture/",
        )
    if "tree" in lower or "phylogen" in lower or "concatenat" in lower:
        return (
            "invalidated",
            "legacy tree inherited faulty inputs and/or fixed model and unsupported 10,000-bootstrap policy",
            f"{base}/results/trees/{run_id}/ and {base}/results/supplement/{run_id}/",
        )
    if "pca" in lower:
        return (
            "invalidated",
            "ordination inherited legacy sample and variant filters",
            f"{base}/results/pca/{run_id}/",
        )
    if "haplotype" in lower:
        return (
            "invalidated",
            "haplotypes inherited legacy callable consensus",
            f"{base}/results/haplotypes/{run_id}/",
        )
    if "reference_verification" in lower or "mtdna_investigation" in lower or "cpdna_investigation" in lower or "analysis_masks" in lower:
        return (
            "invalidated",
            "reference/mask evidence was incomplete, copied, or produced by crashing/permissive code",
            f"{base}/references/",
        )
    if "all_sample_alignment" in lower or "pilot_alignment" in lower or "downstream_sample" in lower:
        return (
            "invalidated",
            "permissive mapping and noncanonical/manual sample-set selection",
            f"{base}/results/qc/{run_id}/ and {base}/metadata/qc/{run_id}/",
        )
    if "review_response" in lower or "submission_figures" in lower or "genome_maps" in lower:
        return (
            "invalidated",
            "figure/report depended on invalidated legacy results",
            f"{base}/reports/",
        )
    return (
        "archived_only",
        "pre-remediation artifact is not a canonical input or deliverable",
        "retired; no direct replacement",
    )


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    run_id = args.run_id
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    validate_publication_config(config)
    report_dir = root / "canonical_publication/reports/manuscript_support" / run_id
    table_dir = root / "canonical_publication/reports/tables" / run_id
    provenance_dir = root / "canonical_publication/provenance"
    manifest_dir = provenance_dir / "manifests"
    invalidation_dir = provenance_dir / "invalidation"
    for directory in (report_dir, table_dir, manifest_dir, invalidation_dir):
        directory.mkdir(parents=True, exist_ok=True)

    archive_manifest = root / "archive_noncanonical/2026-08-17_pre_remediation/manifest.tsv"
    archive_rows = read_tsv(archive_manifest)
    archive_content_digest = validate_inventory(archive_rows, root)
    run_provenance_dir = provenance_dir / "runs" / run_id
    state_path = run_provenance_dir / "reports.json"
    input_state_paths = sorted(
        path for path in run_provenance_dir.rglob("*.json") if path not in {state_path, run_provenance_dir / "ACCEPTANCE.json"}
    )
    declared_inputs = {
        **runtime_provenance(
            root,
            {
                "bcftools": ("bcftools", "--version"),
                "biopython": ("python", "-c", "import Bio; print(Bio.__version__)"),
                "python": ("python", "--version"),
                "ripgrep": ("rg", "--version"),
            },
        ),
        config_path.relative_to(root).as_posix(): sha256_file(config_path),
        "provenance:manifested_pipeline_code": pipeline_code_digest(root),
        **{
            path.relative_to(root).as_posix(): sha256_file(path)
            for path in (
                root / "canonical_publication/README.md",
                root / "canonical_publication/environment.yml",
                root / "canonical_publication/pyproject.toml",
                root / "canonical_publication/validation_environment.yml",
            )
        },
        archive_manifest.relative_to(root).as_posix(): sha256_file(archive_manifest),
        "archive:snapshot_content": archive_content_digest,
        **{path.relative_to(root).as_posix(): sha256_file(path) for path in input_state_paths},
    }
    upstream_fingerprints = {}
    for path in input_state_paths:
        recorded = json.loads(path.read_text()).get("fingerprint")
        digest = recorded.get("digest") if isinstance(recorded, dict) else recorded
        if isinstance(digest, str):
            upstream_fingerprints[path.relative_to(run_provenance_dir).as_posix()] = digest
    fingerprint = build_stage_fingerprint_from_hashes(
        "reports_and_acceptance",
        declared_inputs,
        upstream_fingerprints,
        ["rebuild reports; account for every legacy artifact; verify acceptance; checksum canonical deliverables"],
    )
    if args.resume and state_path.exists():
        saved = json.loads(state_path.read_text())
        validate_resume(saved["fingerprint"]["digest"], fingerprint)
        validate_saved_outputs(root, saved)
        print("resume-valid reports and acceptance")
        return 0
    if state_path.exists():
        raise RuntimeError("Report state already exists; use --resume or a new run ID")
    expected_report_outputs = [
        report_dir / "canonical_analysis_report.md",
        table_dir / "canonical_summary.tsv",
        table_dir / "review_response_method_resolutions.tsv",
        invalidation_dir / f"{run_id}.legacy_artifact_invalidation.tsv",
        manifest_dir / f"{run_id}.final_artifacts.tsv",
        run_provenance_dir / "ACCEPTANCE.json",
        root / "canonical_publication/CURRENT_RUN",
    ]
    if any(path.exists() for path in expected_report_outputs):
        raise RuntimeError("Existing unvalidated report output; preserve it and use a new run ID")

    invalidation_path = invalidation_dir / f"{run_id}.legacy_artifact_invalidation.tsv"
    with invalidation_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "original_path",
                "archived_path",
                "legacy_sha256",
                "disposition",
                "faulty_dependency_or_reason",
                "canonical_replacement",
            ]
        )
        for row in archive_rows:
            disposition, reason, replacement = invalidation_for(row["original_path"], run_id)
            writer.writerow(
                [
                    row["original_path"],
                    row["archived_path"],
                    row["sha256"],
                    disposition,
                    reason,
                    replacement,
                ]
            )

    qc_state = json.loads((root / "canonical_publication/provenance/runs" / run_id / "qc.json").read_text())
    sample_counts = qc_state["sample_counts"]
    high_confidence_variant_counts = {}
    variant_counts = {}
    for organelle in ("chloroplast", "mitochondria"):
        for target, suffix in (
            (high_confidence_variant_counts, "high_confidence_variant_sites"),
            (variant_counts, "primary"),
        ):
            vcf = root / "canonical_publication/results/variants" / run_id / f"{organelle}.{suffix}.vcf.gz"
            target[organelle] = int(
                subprocess.run(
                    ["bcftools", "index", "-n", str(vcf)],
                    capture_output=True,
                    text=True,
                    check=True,
                ).stdout.strip()
            )
    pi_ranges = {}
    population_counts = {}
    pair_counts = {}
    for organelle in ("chloroplast", "mitochondria"):
        populations = read_tsv(root / "canonical_publication/results/popgen" / run_id / f"{organelle}.population_summary.tsv")
        pairs = read_tsv(root / "canonical_publication/results/popgen" / run_id / f"{organelle}.pairwise_hudson_fst.tsv")
        values = [float(row["nucleotide_diversity"]) for row in populations if math.isfinite(float(row["nucleotide_diversity"]))]
        pi_ranges[organelle] = (min(values), max(values)) if values else (math.nan, math.nan)
        population_counts[organelle] = len(populations)
        pair_counts[organelle] = len(pairs)
    conflict_path = root / "canonical_publication/results/supplement" / run_id / "strongly_supported_organelle_conflicts.tsv"
    conflict_count = max(0, len(conflict_path.read_text().splitlines()) - 1)
    concatenated_contributions = {
        row["partition"]: row
        for row in read_tsv(root / "canonical_publication/results/supplement" / run_id / "concatenated_site_contribution.tsv")
    }
    source_state_summary = json.loads((run_provenance_dir / "source_validation.json").read_text())
    declared_missing_sources = source_state_summary.get("declared_missing_provider_entries", [])
    provider_self_reference_warnings = source_state_summary.get("unverifiable_provider_manifest_self_references", [])
    reference_selection = {
        row["organelle"]: row for row in read_tsv(root / "canonical_publication/references/evidence/reference_selection.tsv")
    }
    reference_similarity = {
        row["comparison"]: row for row in read_tsv(root / "canonical_publication/references/evidence/reference_similarity_summary.tsv")
    }
    mitochondria_repeat = read_tsv(root / "canonical_publication/references/evidence/mitochondria_repeat_summary.tsv")[0]
    annotation_projection = {
        row["organelle"]: row
        for row in read_tsv(root / "canonical_publication/references/evidence/annotation_projection/annotation_projection_summary.tsv")
    }
    orientation_checks = {
        row["organelle"]: row
        for row in read_tsv(root / "canonical_publication/references/evidence/reference_orientation_and_boundary_checks.tsv")
    }
    accession_consistency = read_tsv(root / "canonical_publication/references/evidence/external_accession_sequence_consistency.tsv")
    read_backed_reference = {
        row["metric"]: row["value"]
        for row in read_tsv(root / "canonical_publication/references/evidence" / run_id / "read_backed_reference_validation.tsv")
    }
    mt_annotation_overlap = read_tsv(
        root / "canonical_publication/references/evidence" / run_id / "mitochondria_high_confidence_annotation_overlap.tsv"
    )
    mt_features_with_support = sum(row["overlap_status"] in {"partial", "full"} for row in mt_annotation_overlap)
    mt_features_fully_supported = sum(row["overlap_status"] == "full" for row in mt_annotation_overlap)
    reference_identity_medians = {
        organelle: statistics.median(
            float(row["callable_reference_identity"])
            for row in read_tsv(root / "canonical_publication/results/alignments" / run_id / f"{organelle}.callable_summary.tsv")
        )
        for organelle in ("chloroplast", "mitochondria")
    }
    admixture_states = {
        organelle: json.loads((run_provenance_dir / "admixture" / f"{organelle}.json").read_text())
        for organelle in ("chloroplast", "mitochondria")
    }
    sample_manifest_rows = read_tsv(root / "canonical_publication/metadata/samples/samples.tsv")
    preprocessing_rows = read_tsv(root / "canonical_publication/results/qc" / run_id / "read_preprocessing_summary.tsv")
    input_reads = sum(int(row["input_reads"]) for row in preprocessing_rows)
    passing_reads = sum(int(row["passing_reads"]) for row in preprocessing_rows)
    input_bases = sum(int(row["input_bases"]) for row in preprocessing_rows)
    passing_bases = sum(int(row["passing_bases"]) for row in preprocessing_rows)
    weighted_input_q20 = sum(float(row["input_q20_rate"]) * int(row["input_bases"]) for row in preprocessing_rows) / input_bases
    weighted_passing_q20 = sum(float(row["passing_q20_rate"]) * int(row["passing_bases"]) for row in preprocessing_rows) / passing_bases
    adapter_trimmed_reads = sum(int(row["adapter_trimmed_reads"]) for row in preprocessing_rows)
    median_duplication = statistics.median(float(row["duplication_rate"]) for row in preprocessing_rows)
    incomplete_pairs = [row for row in sample_manifest_rows if row["pair_status"] != "complete"]
    unresolved_complete_samples = [
        row["sample_id"] for row in sample_manifest_rows if row["pair_status"] == "complete" and not row["popcode"]
    ]
    report_path = report_dir / "canonical_analysis_report.md"
    report_path.write_text(
        f"# Canonical Dudleya organelle analysis — {run_id}\n\n"
        "## Status\n\n"
        "This report is generated exclusively from the canonical dependency chain. "
        "Pre-remediation files are quarantined and were not analysis inputs.\n\n"
        "## Sample sets\n\n"
        f"- Chloroplast: {sample_counts['chloroplast']} QC-eligible samples.\n"
        f"- Mitochondria: {sample_counts['mitochondria']} QC-eligible samples.\n"
        f"- Shared intersection: {sample_counts['shared']} samples, used only for concatenation.\n\n"
        f"- Incomplete read pairs excluded before analysis: {len(incomplete_pairs)}.\n"
        f"- Provider-manifest entries declared absent from the deposit: {len(declared_missing_sources)}.\n\n"
        f"- Provider manifests with an unauthenticatable self-checksum entry: "
        f"{len(provider_self_reference_warnings)}; these manifest files remain independently SHA-256 inventoried.\n\n"
        f"- Immutable source files passing the SHA-256 inventory: "
        f"{source_state_summary.get('source_inventory_sha256_pass', 0)}/"
        f"{source_state_summary.get('source_inventory_files', 0)}.\n\n"
        "## Read processing and evidence filters\n\n"
        f"Paired reads were adapter-trimmed and filtered with fastp, using Q{int(config['preprocessing']['qualified_quality_phred'])} "
        f"as the qualified-base threshold, rejecting reads with more than "
        f"{int(config['preprocessing']['maximum_unqualified_base_percent'])}% unqualified bases, and requiring "
        f"length ≥{int(config['preprocessing']['minimum_length'])}. BWA-MEM mapping included read groups. "
        f"Evidence with MAPQ <{int(config['mapping']['minimum_mapping_quality'])} was excluded; unmapped, secondary, "
        "supplementary, QC-failed, and duplicate records were removed. Duplicate marking used paired-read fixmate "
        f"metadata, and downstream depth/pileup evidence required base quality ≥{int(config['mapping']['minimum_base_quality'])}.\n\n"
        f"Across {len(preprocessing_rows)} complete pairs, fastp retained {passing_reads:,}/{input_reads:,} reads "
        f"({passing_reads / input_reads:.2%}) and {passing_bases:,}/{input_bases:,} bases "
        f"({passing_bases / input_bases:.2%}). The base-weighted Q20 fraction increased from "
        f"{weighted_input_q20:.2%} to {weighted_passing_q20:.2%}; {adapter_trimmed_reads:,} reads "
        f"({adapter_trimmed_reads / input_reads:.2%}) were adapter-trimmed. The median fastp duplication-rate "
        f"diagnostic was {median_duplication:.2%}; alignment duplicates were subsequently removed before inference.\n\n"
        "Eligibility is organelle-specific and requires at least "
        f"{float(config['qc']['minimum_breadth']) * 100:.0f}% breadth at "
        f"DP≥{int(config['qc']['eligibility_depth'])} over the regenerated organelle unique-mappability mask; "
        "full-reference breadth is reported separately.\n\n"
        "## References and callable masks\n\n"
        f"- External FASTA/GenBank sequence-version consistency: "
        f"{sum(row['sequence_identical'] == 'yes' for row in accession_consistency)}/{len(accession_consistency)} accession pairs PASS.\n"
        f"- Chloroplast: {int(reference_selection['chloroplast']['selected_length']):,} bp after a "
        f"self-BLAST-validated {int(reference_selection['chloroplast']['removed_length']):,} bp redundant terminal-copy trim; "
        f"{float(reference_similarity['chloroplast_selected_vs_NC_085682.1']['query_coverage']):.2%} of selected query bases "
        f"align to NC_085682.1 at {float(reference_similarity['chloroplast_selected_vs_NC_085682.1']['weighted_identity_percent']):.3f}% "
        "position-assigned HSP identity. Both IR copies are excluded only from the mappability denominator; one duplicate "
        "IR copy is excluded from population analysis.\n"
        f"- Mitochondria: the {int(reference_selection['mitochondria']['selected_length']):,} bp candidate is retained intact, but "
        f"only {float(reference_similarity['mitochondria_selected_vs_PV256627.1']['query_coverage']):.2%} of query bases and "
        f"{float(reference_similarity['mitochondria_selected_vs_PV256627.1']['subject_coverage']):.2%} of PV256627.1 are covered "
        f"by qualifying alignments. Self-repeats mask {int(mitochondria_repeat['repeat_masked_bases']):,} bp "
        f"({float(mitochondria_repeat['repeat_masked_fraction']):.2%}); "
        f"{int(read_backed_reference['mitochondria_high_confidence_unique_bases']):,} read-supported unique bases remain in "
        "the final high-confidence mask. Among callable consensus bases, the median eligible-sample identity to the selected "
        f"mitochondrial reference is {reference_identity_medians['mitochondria']:.4%}; this mapping-conditioned concordance is "
        "not an independent assembly validation. The median repeat-mask/unique-site depth ratio is "
        f"{float(read_backed_reference['mitochondria_median_sample_repeat_to_unique_mean_depth_ratio']):.3f}; repeat coordinates "
        "remain excluded regardless of depth. This structural discordance limits whole-mitogenome interpretation.\n"
        "- The dominant external-reference HSP has "
        f"{orientation_checks['chloroplast']['dominant_external_hsp_orientation']} orientation for chloroplast and "
        f"{orientation_checks['mitochondria']['dominant_external_hsp_orientation']} orientation for mitochondria. "
        "After assigning overlapping query positions to the highest-bitscore HSP, the same-orientation fractions among "
        f"covered query bases are {float(orientation_checks['chloroplast']['same_orientation_fraction_of_covered_query']):.2%} "
        f"and {float(orientation_checks['mitochondria']['same_orientation_fraction_of_covered_query']):.2%}, respectively. "
        "These are local alignment diagnostics, not evidence of global collinearity.\n"
        "- External GenBank feature projection recovered "
        f"{annotation_projection['chloroplast']['mapped_features']}/{annotation_projection['chloroplast']['candidate_features']} "
        "chloroplast features and "
        f"{annotation_projection['mitochondria']['mapped_features']}/{annotation_projection['mitochondria']['candidate_features']} "
        "mitochondrial features. Of the projected mitochondrial features, "
        f"{mt_features_with_support}/{len(mt_annotation_overlap)} overlap the final read-backed high-confidence mask and "
        f"{mt_features_fully_supported} are fully contained. These annotations are explicitly draft projections, not de novo annotations.\n"
        f"- The median eligible-sample boundary/interior depth ratio is "
        f"{float(read_backed_reference['chloroplast_median_sample_boundary_to_interior_ratio']):.3f} for chloroplast and "
        f"{float(read_backed_reference['mitochondria_median_sample_boundary_to_interior_ratio']):.3f} for mitochondria; "
        "these values are evidence diagnostics, not proof of circularity.\n\n"
        "## Variants and population statistics\n\n"
        f"Genotypes are haploid and are masked when DP <{int(config['variants']['minimum_depth'])}, "
        f"GQ <{int(config['variants']['minimum_genotype_quality'])}, or either field is missing. Accepted biallelic "
        f"SNP sites require QUAL ≥{int(config['variants']['minimum_site_quality'])} and no more than "
        f"{float(config['variants']['maximum_missing_fraction']):.0%} missing genotypes. Fixed-alternate accepted sites "
        "are retained for consensus; segregating primary summaries include MAC≥1, while PCA and supplementary "
        "ADMIXTURE alone require MAC≥2. The per-input pileup depth cap is "
        f"{int(config['variants']['maximum_per_file_pileup_depth'])} reads per site.\n\n"
        f"- Chloroplast high-confidence variant sites used for consensus (including fixed alternate): "
        f"{high_confidence_variant_counts['chloroplast']}.\n"
        f"- Chloroplast primary variants (including singletons): {variant_counts['chloroplast']}.\n"
        f"- Mitochondrial high-confidence variant sites used for consensus (including fixed alternate): "
        f"{high_confidence_variant_counts['mitochondria']}.\n"
        f"- Mitochondrial primary variants (including singletons): {variant_counts['mitochondria']}.\n"
        f"- Chloroplast callable-site π range: {pi_ranges['chloroplast'][0]:.8g}–{pi_ranges['chloroplast'][1]:.8g}.\n"
        f"- Mitochondrial callable-site π range: {pi_ranges['mitochondria'][0]:.8g}–{pi_ranges['mitochondria'][1]:.8g}.\n"
        "- Pairwise differentiation is signed Hudson ratio-of-sums FST with 1 kb block-bootstrap intervals.\n\n"
        "Supplementary ADMIXTURE tested K=1–12 with ten fixed seeds per K. The minimum mean cross-validation error selected "
        f"K={admixture_states['chloroplast']['selected_k']} for chloroplast "
        f"({'boundary optimum' if admixture_states['chloroplast']['boundary_optimum'] else 'interior optimum'}) and "
        f"K={admixture_states['mitochondria']['selected_k']} for mitochondria "
        f"({'boundary optimum' if admixture_states['mitochondria']['boundary_optimum'] else 'interior optimum'}). "
        "These are sensitivity results under the limitations stated below.\n\n"
        "## Phylogenetic interpretation\n\n"
        "Separate unrooted chloroplast and mitochondrial ModelFinder trees with "
        f"{int(config['phylogeny']['shalrt_replicates']):,} SH-aLRT and "
        f"{int(config['phylogeny']['ultrafast_bootstrap_replicates']):,} ultrafast-bootstrap "
        "replicates are primary. The partitioned concatenated tree is "
        "supplementary. The coordinate-padded chloroplast partition spans "
        f"{float(concatenated_contributions['chloroplast']['coordinate_span_fraction']):.1%} of concatenated coordinates, "
        "but masked all-N coordinates are not treated as evidence: chloroplast contributes "
        f"{float(concatenated_contributions['chloroplast']['at_least_two_callable_site_fraction']):.1%} of sites with at least "
        "two callable shared samples, "
        f"{float(concatenated_contributions['chloroplast']['jointly_callable_site_fraction']):.1%} of jointly callable "
        "shared-sample sites, "
        f"{float(concatenated_contributions['chloroplast']['variable_site_fraction']):.1%} of variable sites, including singletons, "
        "and "
        f"{float(concatenated_contributions['chloroplast']['parsimony_informative_site_fraction']):.1%} "
        "of parsimony-informative sites; "
        "mitochondria contributes "
        f"{float(concatenated_contributions['mitochondria']['at_least_two_callable_site_fraction']):.1%}, "
        f"{float(concatenated_contributions['mitochondria']['jointly_callable_site_fraction']):.1%}, "
        f"{float(concatenated_contributions['mitochondria']['variable_site_fraction']):.1%} and "
        f"{float(concatenated_contributions['mitochondria']['parsimony_informative_site_fraction']):.1%}, respectively. "
        f"{conflict_count} pairs of strongly supported incompatible organelle splits were detected.\n\n"
        "## Limitations\n\n"
        "No nuclear decoy assembly is available. MAPQ, base-quality, depth, and duplicate filters "
        "reduce but cannot eliminate NUMT/NUPT ambiguity. PCA describes linked organelle haplotype variation rather "
        "than independent loci. Diversity, differentiation, and haplotype estimates are conditional on the "
        "organelle-specific cohort and callable masks; mitochondrial estimates do not represent the repeat-rich "
        "whole candidate assembly. ADMIXTURE uses pseudo-diploid linked "
        "organelle markers and is descriptive supplementary clustering, not ancestry inference. "
        "Trees remain unrooted because no defensible outgroup was supplied.\n"
    )
    summary_table = table_dir / "canonical_summary.tsv"
    with summary_table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "organelle",
                "qc_eligible_samples",
                "high_confidence_variant_sites_including_fixed_alternate",
                "primary_variants_including_singletons",
                "population_count",
                "pairwise_fst_count",
                "minimum_callable_pi",
                "maximum_callable_pi",
            ]
        )
        for organelle in ("chloroplast", "mitochondria"):
            writer.writerow(
                [
                    organelle,
                    sample_counts[organelle],
                    high_confidence_variant_counts[organelle],
                    variant_counts[organelle],
                    population_counts[organelle],
                    pair_counts[organelle],
                    f"{pi_ranges[organelle][0]:.12g}",
                    f"{pi_ranges[organelle][1]:.12g}",
                ]
            )
    review_response_table = table_dir / "review_response_method_resolutions.tsv"
    with review_response_table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["legacy_or_review_issue", "canonical_resolution", "evidence"])
        writer.writerows(
            [
                [
                    "forced shared sample count",
                    "independent chloroplast and mitochondrial DP5 breadth eligibility",
                    f"cp={sample_counts['chloroplast']};mt={sample_counts['mitochondria']};shared={sample_counts['shared']}",
                ],
                [
                    "reference-filled uncertain sites and loss of fixed-alternate consensus alleles",
                    "uncertain genotypes are N; consensus accepts every high-confidence SNP site before segregating-site MAC filters",
                    "high-confidence, primary MAC>=1, and ordination MAC>=2 VCF layers",
                ],
                [
                    "SNP-only nucleotide-diversity denominator",
                    "all positions jointly callable within each population",
                    "population summaries plus independent scikit-allel cross-check",
                ],
                [
                    "clamped custom differentiation statistic",
                    "signed Hudson ratio-of-sums with 1 kb block bootstrap",
                    "pairwise numerator, denominator, callable sites, and interval columns",
                ],
                [
                    "10,000-bootstrap requirement",
                    "retired; fixed 1,000 SH-aLRT and 1,000 UFBoot with BNNI",
                    "primary organelle tree states and fixed-seed reproducibility table",
                ],
                [
                    "unsupported ABAB/ABMU rooting alternatives",
                    "retired; primary organelle trees remain unrooted without a defensible outgroup",
                    "tree states, taxon checks, and unrooted publication figures",
                ],
                [
                    "circular mitochondrial genome-map presentation",
                    "retired; the mitochondrial reference is drawn linearly because circularity is not established",
                    "boundary-depth diagnostics and canonical reference/callability map",
                ],
                [
                    "whole mitochondrial candidate treated as equivalent to callable sequence",
                    "self-repeats excluded, followed by an eligible-sample read-backed high-confidence mask",
                    "repeat support, annotation overlap, reference concordance, and mask tables",
                ],
                [
                    "concatenated partition contribution based on padded coordinate span",
                    "report coordinate span separately from callable, variable, and parsimony-informative contributions",
                    "concatenated_site_contribution.tsv",
                ],
                [
                    "primary ADMIXTURE interpretation",
                    "supplementary descriptive sensitivity analysis only",
                    "K=1-12, ten seeds per K, boundary flags, and limitations report",
                ],
            ]
        )

    acceptance_errors = []
    expected_mapping_states = {
        run_provenance_dir / "mapping" / f"{row['sample_id']}.json"
        for row in sample_manifest_rows
        if row["pair_status"] == "complete" and row["analysis_eligible"] == "yes"
    }
    observed_mapping_states = set((run_provenance_dir / "mapping").glob("*.json"))
    if observed_mapping_states != expected_mapping_states:
        acceptance_errors.append("mapping provenance sample IDs disagree with the complete analysis-eligible source manifest")
    required_state_paths = [
        run_provenance_dir / name
        for name in (
            "source_validation.json",
            "references.json",
            "metadata.json",
            "mapping_complete.json",
            "qc.json",
        )
    ]
    required_state_paths.extend(expected_mapping_states)
    for stage in ("variants", "consensus", "pca", "haplotypes", "popgen", "admixture"):
        required_state_paths.extend(run_provenance_dir / stage / f"{organelle}.json" for organelle in ("chloroplast", "mitochondria"))
    required_state_paths.extend(
        [
            run_provenance_dir / "trees/chloroplast.json",
            run_provenance_dir / "trees/mitochondria.json",
            run_provenance_dir / "trees/concatenated.json",
            run_provenance_dir / "trees/conflicts.json",
            run_provenance_dir / "trusted_crosscheck.json",
            run_provenance_dir / "tree_reproducibility.json",
            run_provenance_dir / "figures.json",
        ]
    )
    for required_state_path in required_state_paths:
        if not required_state_path.is_file():
            acceptance_errors.append(f"required stage state is missing: {required_state_path.relative_to(root)}")
            continue
        saved_state = json.loads(required_state_path.read_text())
        if saved_state.get("status") not in {
            "complete",
            "PASS",
            *ACCEPTABLE_SOURCE_VALIDATION_STATUSES,
        }:
            acceptance_errors.append(f"stage did not complete successfully: {required_state_path.relative_to(root)}")
        if "outputs" in saved_state:
            try:
                validate_saved_outputs(root, saved_state)
            except (OSError, KeyError, TypeError, ValueError, RuntimeError) as error:
                acceptance_errors.append(f"stage output validation failed for {required_state_path.relative_to(root)}: {error}")
    if unresolved_complete_samples:
        acceptance_errors.append(f"complete samples lack source-derived population codes: {unresolved_complete_samples[:5]}")
    expected_preprocessing_samples = {
        row["sample_id"] for row in sample_manifest_rows if row["pair_status"] == "complete" and row["analysis_eligible"] == "yes"
    }
    observed_preprocessing_samples = [row["sample_id"] for row in preprocessing_rows]
    if (
        len(observed_preprocessing_samples) != len(set(observed_preprocessing_samples))
        or set(observed_preprocessing_samples) != expected_preprocessing_samples
    ):
        acceptance_errors.append("read-preprocessing summary sample IDs disagree with complete analysis-eligible inputs")
    source_state_path = run_provenance_dir / "source_validation.json"
    if not source_state_path.exists():
        acceptance_errors.append("immutable-source MD5 validation is missing")
    else:
        source_state = json.loads(source_state_path.read_text())
        source_output = root / next(iter(source_state["outputs"]))
        source_validation_rows = read_tsv(source_output)
        if source_state.get("status") not in ACCEPTABLE_SOURCE_VALIDATION_STATUSES or any(
            row["status"] not in {"PASS", "DECLARED_MISSING", "UNVERIFIABLE_SELF_REFERENCE"} for row in source_validation_rows
        ):
            acceptance_errors.append("immutable-source MD5 validation failed")
        if source_state.get("source_inventory_sha256_pass") != source_state.get("source_inventory_files"):
            acceptance_errors.append("immutable-source SHA-256 inventory validation failed")
        provider_status_by_source = {
            row["resolved_source_path"]: row["status"] for row in source_validation_rows if row["resolved_source_path"]
        }
        for sample_row in sample_manifest_rows:
            if sample_row["pair_status"] != "complete" or sample_row["analysis_eligible"] != "yes":
                continue
            for field in ("r1_paths", "r2_paths"):
                for read_path in sample_row[field].split(";"):
                    if provider_status_by_source.get(read_path) != "PASS":
                        acceptance_errors.append(f"analysis-eligible read lacks a passing provider checksum: {read_path}")
    for organelle in ("chloroplast", "mitochondria"):
        expected_pairs = population_counts[organelle] * (population_counts[organelle] - 1) // 2
        if pair_counts[organelle] != expected_pairs:
            acceptance_errors.append(f"{organelle} pair count {pair_counts[organelle]} != {expected_pairs}")
        crosschecks = read_tsv(root / "canonical_publication/results/popgen" / run_id / f"{organelle}.pi_crosscheck.tsv")
        if any(row["exact_match"] != "yes" for row in crosschecks):
            acceptance_errors.append(f"{organelle} pi cross-check failure")
        sample_path = root / "canonical_publication/metadata/qc" / run_id / f"{organelle}_samples.tsv"
        sample_ids = {row["sample_id"] for row in read_tsv(sample_path)}
        alignment_path = root / "canonical_publication/results/alignments" / run_id / f"{organelle}.callable_alignment.fa"
        alignment_ids = {name for name, _ in read_single_fasta(alignment_path, expected_records=len(sample_ids))}
        vcf_sample_sets = []
        for suffix in ("high_confidence_variant_sites", "primary"):
            vcf_path = root / "canonical_publication/results/variants" / run_id / f"{organelle}.{suffix}.vcf.gz"
            vcf_sample_sets.append(
                set(
                    subprocess.run(
                        ["bcftools", "query", "-l", str(vcf_path)],
                        capture_output=True,
                        text=True,
                        check=True,
                    ).stdout.splitlines()
                )
            )
        if sample_ids != alignment_ids or any(sample_ids != vcf_ids for vcf_ids in vcf_sample_sets):
            acceptance_errors.append(f"{organelle} sample IDs disagree")
        tree_path = root / "canonical_publication/results/trees" / run_id / f"{organelle}.primary.treefile"
        tree_ids = [terminal.name for terminal in Phylo.read(tree_path, "newick").get_terminals()]
        if len(tree_ids) != len(set(tree_ids)) or set(tree_ids) != sample_ids:
            acceptance_errors.append(f"{organelle} tree terminal IDs disagree with QC samples")
    shared_sample_ids = {row["sample_id"] for row in read_tsv(root / "canonical_publication/metadata/qc" / run_id / "shared_samples.tsv")}
    concatenated_tree = Phylo.read(
        root / "canonical_publication/results/supplement" / run_id / "concatenated.partitioned.treefile",
        "newick",
    )
    concatenated_tree_ids = [terminal.name for terminal in concatenated_tree.get_terminals()]
    if len(concatenated_tree_ids) != len(set(concatenated_tree_ids)) or set(concatenated_tree_ids) != shared_sample_ids:
        acceptance_errors.append("concatenated tree terminal IDs disagree with shared QC samples")
    trusted_state_path = run_provenance_dir / "trusted_crosscheck.json"
    if not trusted_state_path.exists():
        acceptance_errors.append("trusted scikit-allel cross-check is missing")
    else:
        trusted_state = json.loads(trusted_state_path.read_text())
        trusted_output = root / trusted_state["output"]
        if (
            trusted_state.get("status") != "PASS"
            or sha256_file(trusted_output) != trusted_state.get("output_sha256")
            or any(row["match"] != "yes" for row in read_tsv(trusted_output))
        ):
            acceptance_errors.append("trusted scikit-allel cross-check failed")
    treecheck_path = run_provenance_dir / "tree_reproducibility.json"
    if not treecheck_path.exists():
        acceptance_errors.append("fixed-seed tree reproducibility check is missing")
    else:
        treecheck = json.loads(treecheck_path.read_text())
        treecheck_output = root / next(iter(treecheck["outputs"]))
        if treecheck.get("status") != "PASS" or any(row["strong_topology_reproduced"] != "yes" for row in read_tsv(treecheck_output)):
            acceptance_errors.append("fixed-seed tree reproducibility failed")
    figure_manifest_path = root / "canonical_publication/reports/figures" / run_id / "figure_manifest.tsv"
    if not figure_manifest_path.is_file():
        acceptance_errors.append("canonical figure manifest is missing")
    else:
        figure_rows = read_tsv(figure_manifest_path)
        expected_figure_rows = 12 * 3
        if len(figure_rows) != expected_figure_rows:
            acceptance_errors.append(f"canonical figure manifest has {len(figure_rows)} rows; expected {expected_figure_rows}")
        for row in figure_rows:
            path = root / row["path"]
            if not path.is_file() or sha256_file(path) != row["sha256"]:
                acceptance_errors.append(f"canonical figure manifest checksum failed: {row['path']}")
    if len(archive_rows) != len(read_tsv(invalidation_path)):
        acceptance_errors.append("invalidation report does not account for every archive row")
    snapshot_root = root / "archive_noncanonical/2026-08-17_pre_remediation/snapshot"
    archived_on_disk = {path.relative_to(root).as_posix() for path in snapshot_root.rglob("*") if path.is_file() or path.is_symlink()}
    archived_in_manifest = {row["archived_path"] for row in archive_rows}
    if archived_on_disk != archived_in_manifest:
        acceptance_errors.append("archive manifest paths do not exactly match the legacy snapshot")
    for path in (root / "canonical_publication").rglob("*"):
        if path.is_symlink() and "archive_noncanonical" in str(path.resolve(strict=False)):
            acceptance_errors.append(f"canonical symlink resolves into archive: {path}")
    absolute_hits = subprocess.run(
        [
            "rg",
            "-l",
            "-e",
            "/" + "home/",
            "-e",
            "/" + "Users/",
            "-e",
            "/" + "tmp/",
            "canonical_publication",
            "README.md",
        ],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    ).stdout.splitlines()
    absolute_hits = [
        path
        for path in absolute_hits
        if not path.startswith("canonical_publication/work/")
        and "/__pycache__/" not in path
        and "/.pytest_cache/" not in path
        and "/.ruff_cache/" not in path
    ]
    if absolute_hits:
        acceptance_errors.append(f"absolute workstation paths remain: {absolute_hits[:5]}")

    acceptance_path = provenance_dir / "runs" / run_id / "ACCEPTANCE.json"
    manifest_path = manifest_dir / f"{run_id}.final_artifacts.tsv"
    acceptance_path.write_text(
        json.dumps(
            {
                "run_id": run_id,
                "status": "PASS" if not acceptance_errors else "FAIL",
                "errors": acceptance_errors,
                "sample_counts": sample_counts,
                "high_confidence_variant_counts": high_confidence_variant_counts,
                "variant_counts": variant_counts,
                "population_counts": population_counts,
                "pairwise_counts": pair_counts,
                "strong_conflict_count": conflict_count,
                "final_manifest": manifest_path.relative_to(root).as_posix(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    if acceptance_errors:
        raise RuntimeError(f"Canonical acceptance failed: {acceptance_errors}")
    root_readme = root / "README.md"
    root_readme.write_text(
        "# Dudleya organelle population genomics\n\n"
        "This repository is organized by scientific status:\n\n"
        "- [`source_data/`](source_data/) — immutable raw reads and reference candidates.\n"
        "- [`canonical_publication/`](canonical_publication/) — repaired code and accepted publication outputs.\n"
        "- [`archive_noncanonical/`](archive_noncanonical/) — preserved pre-remediation artifacts "
        "that must not be used for current inference.\n\n"
        f"Current canonical run: `{run_id}` (acceptance PASS). The sole supported entrypoint is "
        "`canonical_publication/run_pipeline.sh`; see `canonical_publication/CURRENT_RUN` and "
        "`canonical_publication/provenance/runs/` for checksummed status.\n"
    )
    current_run_path = root / "canonical_publication/CURRENT_RUN"
    current_run_path.write_text(f"{run_id}\tPASS\n")
    deliverable_roots = [
        root / "canonical_publication/config",
        root / "canonical_publication/pipeline",
        root / "canonical_publication/references/selected",
        root / "canonical_publication/references/annotations",
        root / "canonical_publication/references/evidence/annotation_projection",
        root / "canonical_publication/references/evidence" / run_id,
        root / "canonical_publication/references/masks" / run_id,
        root / "canonical_publication/metadata/samples",
        root / "canonical_publication/metadata/populations",
        root / "canonical_publication/metadata/qc" / run_id,
        *(
            root / "canonical_publication/results" / category / run_id
            for category in (
                "qc",
                "variants",
                "alignments",
                "popgen",
                "pca",
                "haplotypes",
                "trees",
                "supplement",
            )
        ),
        *(
            root / "canonical_publication/reports" / category / run_id
            for category in (
                "figures",
                "tables",
                "manuscript_support",
            )
        ),
        run_provenance_dir,
    ]
    explicit_deliverables = [
        root / "canonical_publication/README.md",
        root / "canonical_publication/environment.yml",
        root / "canonical_publication/pyproject.toml",
        root / "canonical_publication/run_pipeline.sh",
        root / "canonical_publication/validation_environment.yml",
        current_run_path,
        root_readme,
        archive_manifest,
        root / "canonical_publication/provenance/manifests/source_inputs.tsv",
        root / "canonical_publication/provenance/manifests" / f"{run_id}.provider_md5_validation.tsv",
        invalidation_path,
        *(path for path in (root / "canonical_publication/references/masks").iterdir() if path.is_file()),
        *(path for path in (root / "canonical_publication/references/evidence").iterdir() if path.is_file()),
    ]
    excluded_suffixes = (
        ".bam",
        ".bai",
        ".bcf",
        ".ufboot",
        ".ckp.gz",
        ".log",
        ".tbi",
        ".csi",
        ".fai",
        ".amb",
        ".ann",
        ".bwt",
        ".pac",
        ".sa",
        ".ndb",
        ".nhr",
        ".nin",
        ".njs",
        ".nog",
        ".nos",
        ".not",
        ".nsq",
        ".ntf",
        ".nto",
    )
    deliverables = sorted(
        {
            *explicit_deliverables,
            *(
                path
                for base in deliverable_roots
                for path in base.rglob("*")
                if path.is_file()
                and path.name != ".gitkeep"
                and "logs" not in path.parts
                and "fastp" not in path.parts
                and not {"__pycache__", ".pytest_cache", ".ruff_cache"} & set(path.parts)
                and path.suffix != ".pyc"
                and not path.name.endswith(excluded_suffixes)
                and path != manifest_path
            ),
        }
    )
    with manifest_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["canonical_path", "size_bytes", "sha256", "status"])
        for path in deliverables:
            writer.writerow(
                [
                    path.relative_to(root).as_posix(),
                    path.stat().st_size,
                    sha256_file(path),
                    "archive_audit_manifest" if path == archive_manifest else "canonical",
                ]
            )
    for row in read_tsv(manifest_path):
        manifested_path = root / row["canonical_path"]
        if sha256_file(manifested_path) != row["sha256"]:
            raise RuntimeError(f"Final artifact checksum self-verification failed: {row['canonical_path']}")
    report_outputs = [
        report_path,
        summary_table,
        review_response_table,
        invalidation_path,
        acceptance_path,
        manifest_path,
        current_run_path,
        root_readme,
    ]
    state_path.write_text(
        json.dumps(
            {
                "status": "complete",
                "fingerprint": asdict(fingerprint),
                "outputs": {path.relative_to(root).as_posix(): sha256_file(path) for path in report_outputs},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"canonical acceptance PASS; {len(deliverables)} manifested artifacts")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
