"""Prepare the combined organelle reference and pilot sample table.

This is step 2 of the pipeline. It still does not align reads. It checks that
the cpDNA/mtDNA reference is structurally ready for mapping, records whether
required external tools are available, creates indexes only when the tools are
installed, and chooses a small representative pilot sample set.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


EXPECTED_REFERENCE_LENGTHS = {
    "chloroplast": 150274,
    "mitochondria": 243359,
}

REQUIRED_TOOLS = ("bwa", "samtools")
RECOMMENDED_QC_TOOLS = ("fastp", "fastqc", "multiqc", "bcftools")


class ReferenceValidationError(ValueError):
    """Raised when the combined cpDNA/mtDNA reference is not usable."""


@dataclass(frozen=True)
class ReferenceCheck:
    record: str
    observed_length: int | str
    expected_length: int | str
    status: str
    note: str


@dataclass(frozen=True)
class ToolCheck:
    tool: str
    required_for: str
    status: str
    path: str
    note: str


@dataclass(frozen=True)
class IndexCheck:
    index_type: str
    status: str
    files: str
    note: str


def read_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    """Read FASTA record lengths without external dependencies."""

    lengths: dict[str, int] = {}
    current_name = ""
    current_length = 0

    with fasta_path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name:
                    lengths[current_name] = current_length
                current_name = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line)

    if current_name:
        lengths[current_name] = current_length

    return lengths


def validate_reference_records(
    observed: dict[str, int],
    expected: dict[str, int] = EXPECTED_REFERENCE_LENGTHS,
) -> list[ReferenceCheck]:
    """Confirm the combined reference has the expected cpDNA/mtDNA records."""

    observed_names = set(observed)
    expected_names = set(expected)
    missing = expected_names - observed_names
    extra = observed_names - expected_names
    if missing or extra:
        raise ReferenceValidationError(
            "Combined reference records do not match expected organelles. "
            f"Missing: {sorted(missing) or 'none'}; extra: {sorted(extra) or 'none'}"
        )

    checks: list[ReferenceCheck] = []
    for record in sorted(expected):
        observed_length = observed[record]
        expected_length = expected[record]
        status = "PASS" if observed_length == expected_length else "WARN"
        note = (
            "Length matches expected verified reference."
            if status == "PASS"
            else "Length differs from expected verified reference; review before mapping."
        )
        checks.append(
            ReferenceCheck(
                record=record,
                observed_length=observed_length,
                expected_length=expected_length,
                status=status,
                note=note,
            )
        )
    return checks


def check_tools() -> list[ToolCheck]:
    """Record whether the tools needed for the next steps are installed."""

    checks: list[ToolCheck] = []
    for tool in REQUIRED_TOOLS:
        path = shutil.which(tool) or ""
        checks.append(
            ToolCheck(
                tool=tool,
                required_for="reference_indexing_and_mapping",
                status="FOUND" if path else "MISSING",
                path=path,
                note=(
                    "Available on PATH."
                    if path
                    else "Install before reference indexing or read alignment."
                ),
            )
        )

    for tool in RECOMMENDED_QC_TOOLS:
        path = shutil.which(tool) or ""
        checks.append(
            ToolCheck(
                tool=tool,
                required_for="read_qc_or_variant_calling",
                status="FOUND" if path else "MISSING",
                path=path,
                note=(
                    "Available on PATH."
                    if path
                    else "Recommended before full pipeline execution."
                ),
            )
        )
    return checks


def prepare_indexes(reference_path: Path, tool_checks: list[ToolCheck]) -> list[IndexCheck]:
    """Create reference indexes only when the required tools are installed."""

    available = {check.tool for check in tool_checks if check.status == "FOUND"}
    checks: list[IndexCheck] = []

    if "samtools" in available:
        run_command(["samtools", "faidx", reference_path.as_posix()])
        checks.append(
            IndexCheck(
                index_type="samtools_faidx",
                status="CREATED_OR_UPDATED",
                files=f"{reference_path}.fai",
                note="Created with samtools faidx.",
            )
        )
    else:
        checks.append(
            IndexCheck(
                index_type="samtools_faidx",
                status="SKIPPED_TOOL_MISSING",
                files=f"{reference_path}.fai",
                note="samtools is missing; install it before alignment.",
            )
        )

    if "bwa" in available:
        run_command(["bwa", "index", reference_path.as_posix()])
        checks.append(
            IndexCheck(
                index_type="bwa",
                status="CREATED_OR_UPDATED",
                files=";".join(f"{reference_path}{suffix}" for suffix in bwa_suffixes()),
                note="Created with bwa index.",
            )
        )
    else:
        checks.append(
            IndexCheck(
                index_type="bwa",
                status="SKIPPED_TOOL_MISSING",
                files=";".join(f"{reference_path}{suffix}" for suffix in bwa_suffixes()),
                note="bwa is missing; install bwa or adapt the alignment step to bwa-mem2.",
            )
        )

    return checks


def run_command(args: list[str]) -> None:
    subprocess.run(args, check=True)


def bwa_suffixes() -> tuple[str, ...]:
    return (".amb", ".ann", ".bwt", ".pac", ".sa")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_pilot_samples(
    rows: list[dict[str, str]],
    max_samples: int = 15,
) -> list[dict[str, str]]:
    """Choose a deterministic, diverse set of complete paired-end samples."""

    eligible = [
        row
        for row in rows
        if row.get("analysis_status") == "include_primary_paired_end"
        and row.get("r1_paths")
        and row.get("r2_paths")
        and row.get("pair_status") == "complete"
    ]
    selected: list[dict[str, str]] = []
    selected_ids: set[str] = set()

    def add_first(reason: str, predicate) -> None:
        if len(selected) >= max_samples:
            return
        for row in sorted(eligible, key=sample_sort_key):
            sample_id = row["sample_id"]
            if sample_id in selected_ids:
                continue
            if predicate(row):
                selected_ids.add(sample_id)
                selected.append(with_reason(row, reason))
                return

    add_first(
        "representative_initial_du_dash",
        lambda row: row.get("naming_profile") == "initial_du_dash",
    )
    add_first(
        "representative_initial_du_lp",
        lambda row: row.get("naming_profile") == "initial_du_lp",
    )
    add_first("representative_D_cymosa", lambda row: species_group(row) == "D_cymosa")
    add_first(
        "representative_D_abramsii",
        lambda row: species_group(row) == "D_abramsii",
    )
    add_first(
        "representative_D_setchellii",
        lambda row: species_group(row) == "D_setchellii",
    )

    seen_popcodes = {
        row.get("popcode", "") for row in selected if row.get("popcode", "")
    }
    grouped_candidates = main_standard_candidates_by_species(
        eligible, selected_ids, seen_popcodes
    )
    group_order = ("D_cymosa", "D_abramsii", "D_setchellii", "unresolved")
    while len(selected) < max_samples:
        added_in_round = False
        for group in group_order:
            if len(selected) >= max_samples:
                break
            candidates = grouped_candidates.get(group, [])
            while candidates:
                row = candidates.pop(0)
                sample_id = row["sample_id"]
                popcode = row.get("popcode", "")
                if sample_id in selected_ids or popcode in seen_popcodes:
                    continue
                selected_ids.add(sample_id)
                seen_popcodes.add(popcode)
                selected.append(
                    with_reason(row, f"additional_population_{safe_reason(popcode)}")
                )
                added_in_round = True
                break
        if not added_in_round:
            break

    for row in sorted(eligible, key=population_sort_key):
        if len(selected) >= max_samples:
            break
        sample_id = row["sample_id"]
        if sample_id in selected_ids:
            continue
        selected_ids.add(sample_id)
        reason = f"additional_sample_{safe_reason(row.get('popcode') or 'unresolved')}"
        selected.append(with_reason(row, reason))

    return selected


def main_standard_candidates_by_species(
    rows: list[dict[str, str]],
    selected_ids: set[str],
    seen_popcodes: set[str],
) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {
        "D_cymosa": [],
        "D_abramsii": [],
        "D_setchellii": [],
        "unresolved": [],
    }
    for row in sorted(rows, key=population_sort_key):
        sample_id = row["sample_id"]
        popcode = row.get("popcode", "")
        if row.get("naming_profile") != "main_standard" or not popcode:
            continue
        if sample_id in selected_ids or popcode in seen_popcodes:
            continue
        grouped.setdefault(species_group(row), []).append(row)
    return grouped


def with_reason(row: dict[str, str], reason: str) -> dict[str, str]:
    copy = dict(row)
    copy["pilot_reason"] = reason
    return copy


def species_group(row: dict[str, str]) -> str:
    species = row.get("species", "")
    if "cymosa" in species:
        return "D_cymosa"
    if "abramsii" in species:
        return "D_abramsii"
    if "setchellii" in species:
        return "D_setchellii"
    return "unresolved"


def sample_sort_key(row: dict[str, str]) -> tuple[str, str, str]:
    return (
        row.get("batch", ""),
        row.get("popcode", ""),
        row.get("sample_id", ""),
    )


def population_sort_key(row: dict[str, str]) -> tuple[str, str, str, str]:
    return (
        row.get("naming_profile", ""),
        species_group(row),
        row.get("popcode", ""),
        row.get("sample_id", ""),
    )


def safe_reason(value: str) -> str:
    return "".join(character if character.isalnum() else "_" for character in value)


def write_pilot_samples(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else ["sample_id", "pilot_reason"]
    if "pilot_reason" not in fieldnames:
        fieldnames.append("pilot_reason")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_dataclass_tsv(path: Path, rows: Iterable[object], row_type: type[object]) -> None:
    fieldnames = list(row_type.__dataclass_fields__.keys())  # type: ignore[attr-defined]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: getattr(row, field) for field in fieldnames})


def write_summary(
    path: Path,
    reference_path: Path,
    reference_checks: list[ReferenceCheck],
    tool_checks: list[ToolCheck],
    index_checks: list[IndexCheck],
    pilot_rows: list[dict[str, str]],
) -> None:
    missing_tools = [check.tool for check in tool_checks if check.status == "MISSING"]
    lines = [
        "# Step 2 Reference And Pilot Preflight",
        "",
        "This step validates the combined cpDNA/mtDNA reference, records tool",
        "availability, prepares indexes only when tools are installed, and writes",
        "a representative pilot sample table. It does not align reads.",
        "",
        "## Reference",
        "",
        f"- Reference: `{reference_path}`",
        *[
            f"- {check.record}: {check.observed_length} bp "
            f"({check.status}; expected {check.expected_length})"
            for check in reference_checks
        ],
        "",
        "## Tool Availability",
        "",
        *[
            f"- {check.tool}: {check.status}"
            + (f" (`{check.path}`)" if check.path else "")
            for check in tool_checks
        ],
        "",
        "## Index Status",
        "",
        *[
            f"- {check.index_type}: {check.status}. {check.note}"
            for check in index_checks
        ],
        "",
        "## Pilot Sample Set",
        "",
        f"- Pilot samples selected: {len(pilot_rows)}",
        "- Source table: `dudleya_organelle_alignment_pipeline/results/00_manifest/analysis_samples.tsv`",
        "- Missing-mate samples are not eligible for this pilot set.",
        "",
    ]
    if missing_tools:
        lines.extend(
            [
                "## Before Pilot Alignment",
                "",
                "Install or activate an environment containing these missing tools:",
                "",
                *[f"- {tool}" for tool in missing_tools],
                "",
            ]
        )
    path.write_text("\n".join(lines))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Prepare the Dudleya organelle reference and pilot sample table."
    )
    parser.add_argument(
        "--reference",
        type=Path,
        default=Path(
            "dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa"
        ),
        help="Combined cpDNA/mtDNA reference FASTA.",
    )
    parser.add_argument(
        "--analysis-samples",
        type=Path,
        default=Path(
            "dudleya_organelle_alignment_pipeline/results/00_manifest/analysis_samples.tsv"
        ),
        help="Primary paired-end sample table from step 1.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("dudleya_organelle_alignment_pipeline/results/01_reference_pilot"),
        help="Folder for reference checks, tool checks, and pilot samples.",
    )
    parser.add_argument(
        "--max-pilot-samples",
        type=int,
        default=15,
        help="Maximum number of samples to include in the pilot table.",
    )
    parser.add_argument(
        "--skip-indexing",
        action="store_true",
        help="Validate only; do not run samtools faidx or bwa index even if available.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    observed_lengths = read_fasta_lengths(args.reference)
    reference_checks = validate_reference_records(observed_lengths)
    tool_checks = check_tools()
    if args.skip_indexing:
        index_checks = [
            IndexCheck(
                index_type="samtools_faidx",
                status="SKIPPED_BY_USER",
                files=f"{args.reference}.fai",
                note="Indexing skipped by --skip-indexing.",
            ),
            IndexCheck(
                index_type="bwa",
                status="SKIPPED_BY_USER",
                files=";".join(f"{args.reference}{suffix}" for suffix in bwa_suffixes()),
                note="Indexing skipped by --skip-indexing.",
            ),
        ]
    else:
        index_checks = prepare_indexes(args.reference, tool_checks)

    analysis_rows = read_tsv(args.analysis_samples)
    pilot_rows = select_pilot_samples(analysis_rows, args.max_pilot_samples)

    write_dataclass_tsv(
        args.output_dir / "reference_checks.tsv",
        reference_checks,
        ReferenceCheck,
    )
    write_dataclass_tsv(args.output_dir / "tool_checks.tsv", tool_checks, ToolCheck)
    write_dataclass_tsv(args.output_dir / "index_checks.tsv", index_checks, IndexCheck)
    write_pilot_samples(args.output_dir / "pilot_samples.tsv", pilot_rows)
    write_summary(
        args.output_dir / "reference_pilot_summary.md",
        args.reference,
        reference_checks,
        tool_checks,
        index_checks,
        pilot_rows,
    )

    print(f"Reference records checked: {len(reference_checks)}")
    print(f"Pilot samples written: {len(pilot_rows)}")
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
