"""Phase-1 decision records and gates."""

from __future__ import annotations

import json
from pathlib import Path

from .io import read_tsv, write_json, write_tsv


def claim_decision_path(root: Path, run_id: str, *, phase1: bool) -> Path:
    filename = "claim_analysis_decisions.phase1.tsv" if phase1 else "claim_analysis_decisions.tsv"
    return root / f"supplementary_analysis/reports/manuscript_support/{run_id}/{filename}"


def provider_checksum_gate(rows: list[dict[str, str]]) -> tuple[str, list[str]]:
    failures = [row.get("provider_name", "unknown") for row in rows if row["supplementary_status"] == "FAIL"]
    return ("FAIL", failures) if failures else ("PASS", [])


def write_claim_decisions(root: Path, run_id: str) -> list[Path]:
    statuses = read_tsv(root / f"supplementary_analysis/results/sensitivity/{run_id}/sensitivity_status.tsv")
    by_metric: dict[str, list[str]] = {}
    for row in statuses:
        by_metric.setdefault(row["metric"], []).append(row["status"])

    def aggregate(metric: str) -> str:
        values = by_metric.get(metric, [])
        if "FAIL" in values:
            return "FAIL"
        if "PASS_WITH_CAVEAT" in values:
            return "PASS_WITH_CAVEAT"
        return "PASS"

    rows = [
        {
            "claim": "Population diversity estimates are not driven by approved filtering thresholds",
            "analysis": "Filtering and mitochondrial-mask sensitivity",
            "metric": "pi",
            "result_status": aggregate("pi"),
            "required_interpretation_change": "Report ranges and identify threshold-sensitive populations if caveat/fail",
        },
        {
            "claim": "Relative population differentiation is robust to approved filtering thresholds",
            "analysis": "Filtering and mitochondrial-mask sensitivity",
            "metric": "fst",
            "result_status": aggregate("fst"),
            "required_interpretation_change": "Restrict claims to consistently supported pairwise contrasts if caveat/fail",
        },
        {
            "claim": "Leading ordination structure is not a filtering artifact",
            "analysis": "PC1-PC3 Procrustes permutation tests",
            "metric": "pca",
            "result_status": aggregate("pca"),
            "required_interpretation_change": "Describe threshold-dependent structure rather than stable clusters if caveat/fail",
        },
        {
            "claim": "Organelle alignments contain tree-like phylogenetic information",
            "analysis": "Likelihood mapping",
            "metric": "seven_region_likelihood_mapping",
            "result_status": "PENDING_PHASE2",
            "required_interpretation_change": (
                "State insufficient resolution or conflicting signal according to the predeclared decision rule"
            ),
        },
        {
            "claim": "Supported chloroplast and mitochondrial histories can be compared",
            "analysis": "Support-contracted tanglegram and normalized unrooted RF",
            "metric": "supported_topology_compatibility",
            "result_status": "PENDING_PHASE2",
            "required_interpretation_change": "Avoid total-history disagreement language; report numerator and denominator",
        },
        {
            "claim": "Observed population-diversity patterns are not solely consequences of marker count or sample size",
            "analysis": "Site and sample-size resampling",
            "metric": "resampling_distributions",
            "result_status": "PENDING_PHASE2",
            "required_interpretation_change": "Qualify patterns that overlap their predeclared null distributions",
        },
    ]
    output = claim_decision_path(root, run_id, phase1=True)
    write_tsv(output, rows, list(rows[0]), root)
    return [output]


def write_inheritance_evidence(root: Path, run_id: str) -> list[Path]:
    output = root / f"supplementary_analysis/reports/manuscript_support/{run_id}/organelle_inheritance_evidence.md"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "# Organelle inheritance evidence\n\n"
        "No direct experimental study establishing both chloroplast and mitochondrial inheritance in "
        "*Dudleya* or Crassulaceae was identified for this analysis. General angiosperm inheritance "
        "patterns cannot establish the mode in this lineage.\n\n"
        "Permitted manuscript wording:\n\n"
        "> Organelle lineages; inheritance mode was not established in this study.\n\n"
        "Do not substitute ‘maternal lineages,’ ‘seed-mediated lineages,’ or a corresponding dispersal "
        "claim without direct lineage-specific evidence.\n\n"
        "Context reference: Zhang Q, Liu Y, Sodmergen (2003), Examination of the cytoplasmic DNA in "
        "male reproductive cells to determine the potential for cytoplasmic inheritance in 295 "
        "angiosperm species, *Plant and Cell Physiology* 44:941–951, doi:10.1093/pcp/pcg121. "
        "This broad survey is contextual only and does not establish inheritance in Dudleya.\n"
    )
    return [output]


def write_phase1_acceptance(root: Path, run_id: str) -> list[Path]:
    verification = read_tsv(root / f"supplementary_analysis/metadata/qc/{run_id}/metadata_verification.tsv")
    identity = read_tsv(root / f"supplementary_analysis/results/verification/{run_id}/identity/sample_identity_outcomes.tsv")
    claim_path = claim_decision_path(root, run_id, phase1=True)
    inheritance_path = root / f"supplementary_analysis/reports/manuscript_support/{run_id}/organelle_inheritance_evidence.md"
    confirmed = [row["sample_id"] for row in identity if row["outcome"].startswith("confirmed")]
    provider_status, provider_failures = provider_checksum_gate(
        read_tsv(root / f"supplementary_analysis/results/verification/{run_id}/identity/provider_md5_revalidation.tsv")
    )
    metadata_complete = {row["entity"] for row in verification} == {"DUSE", "CY_CAS", "CY_SIE", "TUL2"}
    checks = {
        "metadata_disposition_complete": "PASS" if metadata_complete else "FAIL",
        "resolved_provider_md5_records": provider_status,
        "confirmed_identity_defects_corrected_or_excluded": "PASS" if not confirmed else "FAIL",
        "sensitivity_outputs_present": "PASS"
        if (root / f"supplementary_analysis/results/sensitivity/{run_id}/sensitivity_summary.tsv").is_file()
        else "FAIL",
        "claim_documentation_present": "PASS" if claim_path.is_file() else "FAIL",
        "inheritance_documentation_present": "PASS" if inheritance_path.is_file() else "FAIL",
    }
    payload = {
        "run_id": run_id,
        "status": "PASS" if all(value == "PASS" for value in checks.values()) else "FAIL",
        "checks": checks,
        "confirmed_identity_defects": confirmed,
        "provider_md5_failures": provider_failures,
        "unresolved_identity_samples": sum(row["outcome"] == "unresolved" for row in identity),
    }
    output = root / f"supplementary_analysis/results/verification/{run_id}/phase1_acceptance.json"
    write_json(output, payload, root)
    if payload["status"] != "PASS":
        raise RuntimeError(f"Phase 1 acceptance failed: {json.dumps(checks, sort_keys=True)}")
    return [output]
