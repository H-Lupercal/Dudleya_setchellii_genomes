"""Build the downstream sample set after all-sample organelle alignment QC.

This stage does not call variants or create
alignments. It converts the QC decisions into the sample lists that
variant calling, consensus generation, PCA, tree, Fst, and clustering steps
must use.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


DEFAULT_ANALYSIS_SAMPLE_TABLE = Path(
    "dudleya_organelle_alignment_pipeline/results/00_manifest/analysis_samples.tsv"
)
DEFAULT_UPSTREAM_EXCLUDED_TABLE = Path(
    "dudleya_organelle_alignment_pipeline/results/00_manifest/excluded_samples.tsv"
)
DEFAULT_QC_DECISION_TABLE = Path(
    "dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/"
    "downstream_sample_qc_decisions.tsv"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set"
)
DEFAULT_EXPECTED_INCLUDED = 275


INCLUDED_FIELDNAMES = [
    "sample_id",
    "batch",
    "naming_profile",
    "species",
    "popcode",
    "population_name",
    "du_id",
    "lp_id",
    "r1_paths",
    "r2_paths",
    "downstream_cpDNA_use",
    "downstream_mtDNA_use",
    "include_reason",
]


EXCLUDED_FIELDNAMES = [
    "sample_id",
    "batch",
    "naming_profile",
    "species",
    "popcode",
    "population_name",
    "du_id",
    "lp_id",
    "exclusion_stage",
    "downstream_cpDNA_use",
    "downstream_mtDNA_use",
    "exclusion_reason",
    "evidence",
]


class DownstreamSampleSetError(RuntimeError):
    """Raised when this stage cannot build a safe downstream sample set."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build_downstream_sample_set(
    analysis_rows: list[dict[str, str]],
    upstream_excluded_rows: list[dict[str, str]],
    qc_decision_rows: list[dict[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Split samples into downstream included and excluded rows."""

    qc_decisions = {
        row["sample_id"]: row
        for row in qc_decision_rows
        if row.get("ignored_downstream", "").lower() == "yes"
    }
    analysis_sample_ids = {row["sample_id"] for row in analysis_rows}
    unknown_qc_ids = sorted(set(qc_decisions) - analysis_sample_ids)
    if unknown_qc_ids:
        raise DownstreamSampleSetError(
            "QC decision table contains sample IDs absent from analysis samples: "
            + ", ".join(unknown_qc_ids)
        )

    included: list[dict[str, str]] = []
    excluded: list[dict[str, str]] = []
    for row in analysis_rows:
        sample_id = row["sample_id"]
        if row.get("analysis_status") != "include_primary_paired_end":
            continue
        if sample_id in qc_decisions:
            excluded.append(build_qc_excluded_row(row, qc_decisions[sample_id]))
        else:
            included.append(build_included_row(row))

    for row in upstream_excluded_rows:
        excluded.append(build_upstream_excluded_row(row))

    return included, excluded


def build_included_row(row: dict[str, str]) -> dict[str, str]:
    return {
        "sample_id": row.get("sample_id", ""),
        "batch": row.get("batch", ""),
        "naming_profile": row.get("naming_profile", ""),
        "species": row.get("species", ""),
        "popcode": row.get("popcode", ""),
        "population_name": row.get("population_name", ""),
        "du_id": row.get("du_id", ""),
        "lp_id": row.get("lp_id", ""),
        "r1_paths": row.get("r1_paths", ""),
        "r2_paths": row.get("r2_paths", ""),
        "downstream_cpDNA_use": "include",
        "downstream_mtDNA_use": "include",
        "include_reason": "passes all-sample-alignment downstream QC",
    }


def build_qc_excluded_row(
    analysis_row: dict[str, str],
    qc_decision_row: dict[str, str],
) -> dict[str, str]:
    return {
        "sample_id": analysis_row.get("sample_id", ""),
        "batch": analysis_row.get("batch", ""),
        "naming_profile": analysis_row.get("naming_profile", ""),
        "species": analysis_row.get("species", ""),
        "popcode": analysis_row.get("popcode", ""),
        "population_name": analysis_row.get("population_name", ""),
        "du_id": analysis_row.get("du_id", ""),
        "lp_id": analysis_row.get("lp_id", ""),
        "exclusion_stage": "step5_downstream_qc",
        "downstream_cpDNA_use": qc_decision_row.get("downstream_cpDNA_use", "exclude"),
        "downstream_mtDNA_use": qc_decision_row.get("downstream_mtDNA_use", "exclude"),
        "exclusion_reason": qc_decision_row.get("reason", ""),
        "evidence": qc_decision_row.get("evidence", ""),
    }


def build_upstream_excluded_row(row: dict[str, str]) -> dict[str, str]:
    reason = row.get("pair_status") or row.get("analysis_status") or "upstream_excluded"
    evidence = row.get("analysis_note", "")
    if row.get("r1_count") or row.get("r2_count"):
        evidence = append_evidence(
            evidence,
            f"r1_count={row.get('r1_count', '')};r2_count={row.get('r2_count', '')}",
        )
    return {
        "sample_id": row.get("sample_id", ""),
        "batch": row.get("batch", ""),
        "naming_profile": row.get("naming_profile", ""),
        "species": row.get("species", ""),
        "popcode": row.get("popcode", ""),
        "population_name": row.get("population_name", ""),
        "du_id": row.get("du_id", ""),
        "lp_id": row.get("lp_id", ""),
        "exclusion_stage": "step0_manifest",
        "downstream_cpDNA_use": "exclude",
        "downstream_mtDNA_use": "exclude",
        "exclusion_reason": reason,
        "evidence": evidence,
    }


def append_evidence(base: str, addition: str) -> str:
    if base and addition:
        return f"{base};{addition}"
    return base or addition


def write_downstream_sample_set_outputs(
    included_rows: list[dict[str, str]],
    excluded_rows: list[dict[str, str]],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(output_dir / "included_samples.tsv", included_rows, INCLUDED_FIELDNAMES)
    write_tsv(output_dir / "excluded_samples.tsv", excluded_rows, EXCLUDED_FIELDNAMES)
    write_report(
        output_dir / "downstream_sample_set_report.md",
        included_rows=included_rows,
        excluded_rows=excluded_rows,
    )


def write_report(
    path: Path,
    included_rows: list[dict[str, str]],
    excluded_rows: list[dict[str, str]],
) -> None:
    stage_counts: dict[str, int] = {}
    for row in excluded_rows:
        stage = row.get("exclusion_stage", "unknown")
        stage_counts[stage] = stage_counts.get(stage, 0) + 1

    lines = [
        "# Downstream Sample Set",
        "",
        "This step defines the sample set for downstream haploid variant calling,",
        "consensus FASTA generation, cpDNA/mtDNA all-sample alignments, PCA,",
        "phylogenetic trees, Fst, and structure/admixture-style clustering.",
        "",
        "## Summary",
        "",
        f"- Included samples: {len(included_rows)}",
        f"- Excluded samples: {len(excluded_rows)}",
        "",
        "## Exclusions By Stage",
        "",
    ]
    for stage in sorted(stage_counts):
        lines.append(f"- {stage}: {stage_counts[stage]}")

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- `included_samples.tsv`: samples to use in primary downstream analyses.",
            "- `excluded_samples.tsv`: samples excluded before downstream analyses.",
            "",
            "The included sample set should be used by variant calling and",
            "all later population-genetic outputs.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def generate_downstream_sample_set(
    analysis_sample_table: Path = DEFAULT_ANALYSIS_SAMPLE_TABLE,
    upstream_excluded_table: Path = DEFAULT_UPSTREAM_EXCLUDED_TABLE,
    qc_decision_table: Path = DEFAULT_QC_DECISION_TABLE,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    expected_included: int | None = DEFAULT_EXPECTED_INCLUDED,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    analysis_rows = read_tsv(analysis_sample_table)
    upstream_excluded_rows = read_tsv(upstream_excluded_table)
    qc_decision_rows = read_tsv(qc_decision_table)
    included_rows, excluded_rows = build_downstream_sample_set(
        analysis_rows=analysis_rows,
        upstream_excluded_rows=upstream_excluded_rows,
        qc_decision_rows=qc_decision_rows,
    )
    if expected_included is not None and len(included_rows) != expected_included:
        raise DownstreamSampleSetError(
            f"Expected {expected_included} downstream included samples, "
            f"found {len(included_rows)}."
        )
    write_downstream_sample_set_outputs(included_rows, excluded_rows, output_dir)
    return included_rows, excluded_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build downstream sample include/exclude tables."
    )
    parser.add_argument("--analysis-sample-table", type=Path, default=DEFAULT_ANALYSIS_SAMPLE_TABLE)
    parser.add_argument("--upstream-excluded-table", type=Path, default=DEFAULT_UPSTREAM_EXCLUDED_TABLE)
    parser.add_argument("--qc-decision-table", type=Path, default=DEFAULT_QC_DECISION_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--expected-included",
        type=int,
        default=DEFAULT_EXPECTED_INCLUDED,
        help="Expected number of downstream included samples. Use -1 to disable.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    expected_included = None if args.expected_included < 0 else args.expected_included
    included_rows, excluded_rows = generate_downstream_sample_set(
        analysis_sample_table=args.analysis_sample_table,
        upstream_excluded_table=args.upstream_excluded_table,
        qc_decision_table=args.qc_decision_table,
        output_dir=args.output_dir,
        expected_included=expected_included,
    )
    print(f"Downstream included samples: {len(included_rows)}")
    print(f"Downstream excluded samples: {len(excluded_rows)}")
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
