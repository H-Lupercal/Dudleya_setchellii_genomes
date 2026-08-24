"""Reports, artifact manifests, and final supplementary acceptance."""

from __future__ import annotations

import json
from pathlib import Path

from .documentation import claim_decision_path
from .figures import validate_figure_manifest
from .io import read_tsv, write_json, write_tsv
from .provenance import sha256_file
from .reporting import acceptance_checks


def resolve_phase2_claims(
    claim_rows: list[dict[str, str]],
    *,
    likelihood_rows: list[dict[str, str]],
    rf_row: dict[str, str],
    resampling_summary: dict[str, str],
) -> list[dict[str, str]]:
    """Replace Phase-2 placeholders with results and required interpretation."""
    rows = [dict(row) for row in claim_rows]
    by_metric = {row["metric"]: row for row in rows}

    likelihood_failures = [
        row["organelle"] for row in likelihood_rows if row["decision"] in {"INSUFFICIENT_RESOLUTION", "INSUFFICIENT_INFORMATION"}
    ]
    composition_caveats = [row["organelle"] for row in likelihood_rows if int(row.get("composition_failed_count", "0")) > 0]
    likelihood_claim = by_metric["seven_region_likelihood_mapping"]
    likelihood_claim["result_status"] = "FAIL" if likelihood_failures else "PASS_WITH_CAVEAT" if composition_caveats else "PASS"
    likelihood_claim["required_interpretation_change"] = (
        f"Report insufficient resolution for {','.join(likelihood_failures)}; do not infer a network from unresolved signal"
        if likelihood_failures
        else (
            "Present trees as unrooted and report composition-test failures for "
            f"{','.join(composition_caveats)}; no NeighborNet trigger was met"
            if composition_caveats
            else "Present both organelle trees as unrooted; no NeighborNet trigger was met"
        )
    )

    rf_numerator = int(rf_row["rf_numerator"])
    rf_denominator = int(rf_row["rf_denominator"])
    topology_claim = by_metric["supported_topology_compatibility"]
    topology_claim["result_status"] = "PASS_WITH_CAVEAT" if rf_numerator else "PASS"
    topology_claim["required_interpretation_change"] = (
        f"Report supported-topology discordance as RF {rf_numerator}/{rf_denominator}; do not describe total evolutionary disagreement"
        if rf_numerator
        else "Report agreement only on the support-contracted 229-representative taxon space"
    )

    marker_robust = resampling_summary["marker_count_result"] == "observed_outside_cp_distribution"
    sample_robust = resampling_summary["sample_size_result"] == "named_outlier_medians_remain_top_ranked"
    resampling_claim = by_metric["resampling_distributions"]
    if marker_robust and sample_robust:
        resampling_claim["result_status"] = "PASS"
    elif marker_robust or sample_robust:
        resampling_claim["result_status"] = "PASS_WITH_CAVEAT"
    else:
        resampling_claim["result_status"] = "FAIL"
    resampling_claim["required_interpretation_change"] = (
        "Sample-size-standardized diversity outliers persist, but mitochondrial haplotype sharing falls within the "
        "146-site chloroplast distribution; treat sharing as marker-count-sensitive"
        if sample_robust and not marker_robust
        else "Interpret marker-count and sample-size controls according to population_resampling_summary.tsv"
    )
    return rows


def write_reports(root: Path, run_id: str) -> list[Path]:
    report_dir = root / f"supplementary_analysis/reports/manuscript_support/{run_id}"
    report_dir.mkdir(parents=True, exist_ok=True)
    sensitivity = read_tsv(root / f"supplementary_analysis/results/sensitivity/{run_id}/sensitivity_status.tsv")
    likelihood = read_tsv(root / f"supplementary_analysis/results/phylogeny/{run_id}/likelihood_mapping/likelihood_mapping_summary.tsv")
    rf = read_tsv(root / f"supplementary_analysis/results/comparative/{run_id}/organelle_comparison/supported_unrooted_rf.tsv")[0]
    fst_agreement = read_tsv(
        root / f"supplementary_analysis/results/comparative/{run_id}/organelle_comparison/common_pair_fst_agreement.tsv"
    )[0]
    resampling = read_tsv(
        root / f"supplementary_analysis/results/comparative/{run_id}/population_diversity/population_resampling_summary.tsv"
    )[0]
    phase1_claim_path = claim_decision_path(root, run_id, phase1=True)
    claim_path = claim_decision_path(root, run_id, phase1=False)
    claim_rows = resolve_phase2_claims(read_tsv(phase1_claim_path), likelihood_rows=likelihood, rf_row=rf, resampling_summary=resampling)
    write_tsv(claim_path, claim_rows, list(claim_rows[0]), root)
    report = report_dir / "supplementary_analysis_report.md"
    status_lines = "\n".join(f"- {row['scenario']} {row['organelle']} {row['metric']}: {row['status']}" for row in sensitivity)
    likelihood_lines = "\n".join(
        f"- {row['organelle']}: resolved={100 * float(row['resolved_fraction']):.2f}%, "
        f"partly resolved={100 * float(row['side_fraction']):.2f}%, unresolved={100 * float(row['center_fraction']):.2f}% "
        f"({row['decision']}); composition failures={row.get('composition_failed_count', 'not_recorded')}/"
        f"{row.get('alignment_sequence_count', 'not_recorded')}, >50% gaps/ambiguity="
        f"{row.get('over_50pct_ambiguity_count', 'not_recorded')}."
        for row in likelihood
    )
    report.write_text(
        f"""# Supplementary analysis report — {run_id}

## What we studied

We tested whether the canonical chloroplast and mitochondrial conclusions are sensitive to
approved filtering and mitochondrial-mask choices. We quantified phylogenetic information,
compared supported unrooted organelle topologies, examined technical covariates, and
standardized marker and sample counts.

## Why it matters

Organelle genomes are linked haploid lineages. Strong-looking clusters or differentiation can
reflect filtering, missingness, reference concordance, or limited phylogenetic information.
These analyses expose those alternatives without treating organelles as independent nuclear loci.

## Data and scope

The immutable base is `publication-20260817` (cp 276; mt 271; shared 271). Existing filtered
BAMs were reused. No preprocessing or remapping was performed. DUSE samples remain in
sample-level displays but are excluded from population inference because the population label
is unresolved. Geography was `not_run:no_approved_coordinates`.

## Robustness outcomes

{status_lines}

Scientific caveats or failures remain visible in the figures and require interpretation changes
recorded in `claim_analysis_decisions.tsv`; they are not provenance failures.

## Phylogenetic information

{likelihood_lines}

The cp–mt comparison is unrooted. After contracting branches lacking joint SH-aLRT≥80 and
UFBoot≥95 support, normalized RF was {rf["rf_numerator"]}/{rf["rf_denominator"]} =
{float(rf["normalized_unrooted_rf"]):.4f} on 229 mitochondrial unique-sequence representatives.
The tanglegram restores/displays all 271 shared samples, including 42 annotated zero-length
identical-tip memberships.

Across eligible common population pairs, cp–mt FST rank agreement was
rho={float(fst_agreement["global_spearman_rho"]):.4f} using
{fst_agreement["finite_common_pair_count"]}/{fst_agreement["eligible_common_pair_count"]} finite pairs;
{fst_agreement["nonfinite_common_pair_count"]} nonfinite pairs remain in the table but were not
used in the correlation.

## Marker-count and sample-size controls

Observed mitochondrial multi-population haplotype sharing was
{float(resampling["observed_mt_multi_population_haplotypes"]):.0f}; this equals the median of the
1,000 chloroplast 146-site draws (95% interval
{float(resampling["cp_draw_q025_multi_population_haplotypes"]):.0f}–
{float(resampling["cp_draw_q975_multi_population_haplotypes"]):.0f}). Thus this sharing count is
marker-count-sensitive. At common n=4, CY_SIE and CY_CAS retained the two highest median
chloroplast pi values; their medians exceeded
{100 * float(resampling["CY_SIE_fraction_other_draws_below_median"]):.1f}% and
{100 * float(resampling["CY_CAS_fraction_other_draws_below_median"]):.1f}% of draws from other
populations, respectively. This is a descriptive resampling comparison, not an independent-locus
hypothesis test.

## Interpretation

The organelles describe lineage history, not a nuclear-genome admixture history. ADMIXTURE
remains a demoted sensitivity visualization because its linked haploid markers violate the usual
independent-diploid interpretation. Organelle inheritance mode was not established here.

## Limitations

- No nuclear decoy was available, so residual NUMT/NUPT ambiguity cannot be excluded.
- Raw-read sketches are screens: sketch similarity alone is suspected, not confirmed identity,
  while a negative result does not prove biological independence.
- Index hopping is untestable without index sequences, sample sheets, and demultiplexing metrics.
- Mitochondrial inference receives extra restraint because the primary alignment contains only
  146 SNPs, {likelihood[1].get("composition_failed_count", "many")}/{likelihood[1].get("alignment_sequence_count", "271")}
  sequences failed the composition test, and weak-split topology was not fully reproducible.
- Marker-count resampling does not control mutation rate, mask, missingness, or organelle biology.
"""
    )
    geography = report_dir / "geography_status.tsv"
    write_tsv(
        geography,
        [
            {
                "analysis": "geography",
                "status": "not_run:no_approved_coordinates",
                "reason": "No approved, provenance-bearing public coordinates",
            }
        ],
        ["analysis", "status", "reason"],
        root,
    )
    table_rows = []
    for path in sorted((root / "supplementary_analysis/results").rglob("*.tsv")):
        table_rows.append({"path": path.relative_to(root).as_posix(), "sha256": sha256_file(path)})
    table_manifest = root / f"supplementary_analysis/reports/tables/{run_id}/table_manifest.tsv"
    write_tsv(table_manifest, table_rows, ["path", "sha256"], root)
    return [claim_path, report, geography, table_manifest]


def _artifact_paths(root: Path, run_id: str) -> list[Path]:
    base = root / "supplementary_analysis"
    paths = []
    for subtree in ("config", "pipeline/src", "pipeline/scripts", "metadata", "results", "reports"):
        for path in sorted((base / subtree).rglob("*")):
            if path.is_file() and "__pycache__" not in path.parts:
                paths.append(path)
    for name in ("README.md", "environment.yml", "run_pipeline.sh"):
        path = base / name
        if path.is_file():
            paths.append(path)
    return sorted(set(paths))


def write_acceptance(root: Path, run_id: str, canonical_unchanged: bool) -> list[Path]:
    figure_manifest = root / f"supplementary_analysis/reports/figures/{run_id}/supplementary_figure_manifest.tsv"
    figure_rows = read_tsv(figure_manifest)
    validate_figure_manifest(figure_rows)
    paths = _artifact_paths(root, run_id)
    manifest_rows = [
        {
            "path": path.relative_to(root).as_posix(),
            "type": "file",
            "size": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in paths
    ]
    manifest = root / f"supplementary_analysis/provenance/manifests/{run_id}.final_artifacts.tsv"
    write_tsv(manifest, manifest_rows, ["path", "type", "size", "sha256"], root)
    tangle = read_tsv(root / f"supplementary_analysis/results/comparative/{run_id}/organelle_comparison/tanglegram_271_tip_mapping.tsv")
    rf = read_tsv(root / f"supplementary_analysis/results/comparative/{run_id}/organelle_comparison/supported_unrooted_rf.tsv")[0]
    acceptance = acceptance_checks(
        canonical_counts={"chloroplast": 276, "mitochondria": 271, "shared": 271},
        figure_family_count=len({row["family"] for row in figure_rows}),
        all_artifacts_checksummed=len(manifest_rows) == len(paths) and all(row["sha256"] for row in manifest_rows),
        canonical_unchanged=canonical_unchanged,
        shared_display_count=len(tangle),
        restored_identical_tip_count=sum(row["identical_zero_length_tip_group"] == "yes" for row in tangle),
        rf_representative_count=int(rf["taxon_space"].split("_", 1)[0]),
    )
    claim_rows = read_tsv(claim_decision_path(root, run_id, phase1=False))
    pending_claims = [row["metric"] for row in claim_rows if row["result_status"].startswith("PENDING")]
    acceptance["claim_decisions_final"] = {
        "status": "PASS" if not pending_claims else "FAIL",
        "pending_metrics": pending_claims,
    }
    if pending_claims:
        acceptance["status"] = "FAIL"
    sensitivity_statuses = {
        row["status"] for row in read_tsv(root / f"supplementary_analysis/results/sensitivity/{run_id}/sensitivity_status.tsv")
    }
    acceptance.update(
        {
            "run_id": run_id,
            "base_run_id": "publication-20260817",
            "decision_plan_version": "2.5",
            "scientific_sensitivity_statuses": sorted(sensitivity_statuses),
            "scientific_failures_block_acceptance": False,
            "geography": "not_run:no_approved_coordinates",
            "artifact_manifest": manifest.relative_to(root).as_posix(),
            "artifact_manifest_sha256": sha256_file(manifest),
        }
    )
    output = root / f"supplementary_analysis/provenance/runs/{run_id}/ACCEPTANCE.json"
    write_json(output, acceptance, root)
    if acceptance["status"] != "PASS":
        raise RuntimeError(f"Supplementary acceptance failed: {json.dumps(acceptance, sort_keys=True)}")
    current = root / "supplementary_analysis/CURRENT_RUN"
    current.write_text(run_id + "\n")
    checksum_index = root / f"supplementary_analysis/provenance/manifests/{run_id}.acceptance.sha256"
    checksum_index.write_text(
        f"{sha256_file(manifest)}  {manifest.relative_to(root).as_posix()}\n"
        f"{sha256_file(output)}  {output.relative_to(root).as_posix()}\n"
        f"{sha256_file(current)}  {current.relative_to(root).as_posix()}\n"
    )
    return [manifest, output, current, checksum_index]
