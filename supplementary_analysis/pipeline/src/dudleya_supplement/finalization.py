"""Reports, artifact manifests, and final supplementary acceptance."""

from __future__ import annotations

import json
from pathlib import Path

from .figures import validate_figure_manifest
from .io import read_tsv, write_json, write_tsv
from .provenance import sha256_file
from .reporting import acceptance_checks


def write_reports(root: Path, run_id: str) -> list[Path]:
    report_dir = root / f"supplementary_analysis/reports/manuscript_support/{run_id}"
    report_dir.mkdir(parents=True, exist_ok=True)
    sensitivity = read_tsv(root / f"supplementary_analysis/results/sensitivity/{run_id}/sensitivity_status.tsv")
    likelihood = read_tsv(root / f"supplementary_analysis/results/phylogeny/{run_id}/likelihood_mapping/likelihood_mapping_summary.tsv")
    rf = read_tsv(root / f"supplementary_analysis/results/comparative/{run_id}/organelle_comparison/supported_unrooted_rf.tsv")[0]
    report = report_dir / "supplementary_analysis_report.md"
    status_lines = "\n".join(f"- {row['scenario']} {row['organelle']} {row['metric']}: {row['status']}" for row in sensitivity)
    likelihood_lines = "\n".join(
        f"- {row['organelle']}: resolved={100 * float(row['resolved_fraction']):.2f}%, "
        f"partly resolved={100 * float(row['side_fraction']):.2f}%, unresolved={100 * float(row['center_fraction']):.2f}% "
        f"({row['decision']})."
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
  146 SNPs and weak-split topology was not fully reproducible.
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
    return [report, geography, table_manifest]


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
        rf_representative_count=int(rf["taxon_space"].split("_", 1)[0]),
    )
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
