"""Build the first pipeline input manifest for Dudleya organelle alignments.

This step does not align reads. It answers the boring but critical questions
that make alignment safe: which FASTQ files exist, which R1/R2 files pair
together, which naming convention each sample follows, and which samples have
population metadata.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


FASTQ_NAME_RE = re.compile(
    r"^(?P<prefix>.+)_"
    r"(?P<sequencing_sample>S\d+)_"
    r"(?P<lane>L\d+)_"
    r"(?P<read>R[12])_"
    r"(?P<chunk>\d+)"
    r"(?P<extension>\.(?:fastq|fq)(?:-\d+)?(?:\.gz)?)$",
    re.IGNORECASE,
)

FASTQ_FILE_RE = re.compile(r"\.(?:fastq|fq)(?:-\d+)?(?:\.gz)?$", re.IGNORECASE)

MAIN_STANDARD_RE = re.compile(
    r"^(?P<popcode>.+)_"
    r"(?P<lp_id>LP_\d+)[_-]"
    r"(?P<du_id>Du-[A-Za-z0-9]+)$"
)

INITIAL_DU_LP_RE = re.compile(r"^(?P<du_id>DU\d+)(?P<lp_id>LP\d+)$")
INITIAL_DU_DASH_RE = re.compile(r"^(?P<du_id>DU-[A-Za-z0-9]+)$")


@dataclass(frozen=True)
class PopulationCode:
    code: str
    species: str
    population_name: str


@dataclass(frozen=True)
class FastqRecord:
    path: Path
    filename: str
    batch: str
    sample_id: str
    naming_profile: str
    sequencing_sample: str
    lane: str
    read: str
    chunk: str
    popcode: str = ""
    du_id: str = ""
    lp_id: str = ""


@dataclass(frozen=True)
class ManifestRow:
    sample_id: str
    batch: str
    naming_profile: str
    popcode: str
    species: str
    population_name: str
    du_id: str
    lp_id: str
    sequencing_samples: str
    lanes: str
    r1_paths: str
    r2_paths: str
    r1_count: int
    r2_count: int
    pair_status: str
    metadata_status: str
    analysis_status: str
    analysis_note: str


@dataclass(frozen=True)
class ManifestIssue:
    sample_id: str
    batch: str
    issue_type: str
    details: str


def clean_cell(value: str | None) -> str:
    """Normalize spreadsheet cells without changing meaningful IDs."""

    if value is None:
        return ""
    return value.strip().strip('"').strip()


def discover_fastq_paths(data_root: Path) -> list[Path]:
    """Find FASTQ files under the data root.

    This uses Python filesystem traversal, not `git`, because the raw data
    folder is intentionally ignored by `.gitignore`.
    """

    paths: list[Path] = []
    for path in data_root.rglob("*"):
        if not path.is_file():
            continue
        if FASTQ_FILE_RE.search(path.name):
            paths.append(path)
    return sorted(paths)


def infer_batch(path: Path) -> str:
    """Return the sequencing-date/results-folder pair when present."""

    parts = path.parts
    if "genomicsDrive_data_dump" in parts:
        start = parts.index("genomicsDrive_data_dump") + 1
        if len(parts) > start + 1:
            return f"{parts[start]}/{parts[start + 1]}"
    if len(parts) >= 3:
        return f"{parts[-3]}/{parts[-2]}"
    return path.parent.as_posix()


def classify_prefix(prefix: str) -> tuple[str, str, str, str]:
    """Classify a filename prefix into known Dudleya naming profiles.

    Returns naming_profile, popcode, du_id, and lp_id. The initial genome
    batches do not encode population metadata in the filename, so their
    popcode is intentionally blank.
    """

    main_match = MAIN_STANDARD_RE.match(prefix)
    if main_match:
        return (
            "main_standard",
            main_match.group("popcode"),
            main_match.group("du_id"),
            main_match.group("lp_id"),
        )

    du_lp_match = INITIAL_DU_LP_RE.match(prefix)
    if du_lp_match:
        return (
            "initial_du_lp",
            "",
            du_lp_match.group("du_id"),
            du_lp_match.group("lp_id"),
        )

    du_dash_match = INITIAL_DU_DASH_RE.match(prefix)
    if du_dash_match:
        return ("initial_du_dash", "", du_dash_match.group("du_id"), "")

    return ("unrecognized", "", "", "")


def parse_fastq_path(path: Path) -> FastqRecord:
    """Parse one Illumina-style FASTQ filename into pipeline metadata."""

    match = FASTQ_NAME_RE.match(path.name)
    if not match:
        raise ValueError(f"Not a recognized Illumina FASTQ name: {path}")

    prefix = match.group("prefix")
    naming_profile, popcode, du_id, lp_id = classify_prefix(prefix)

    return FastqRecord(
        path=path,
        filename=path.name,
        batch=infer_batch(path),
        sample_id=prefix,
        naming_profile=naming_profile,
        sequencing_sample=match.group("sequencing_sample"),
        lane=match.group("lane"),
        read=match.group("read").upper(),
        chunk=match.group("chunk"),
        popcode=popcode,
        du_id=du_id,
        lp_id=lp_id,
    )


def load_population_codes(csv_path: Path) -> dict[str, PopulationCode]:
    """Load Evan's population-code spreadsheet as a code-indexed dictionary."""

    if not csv_path.exists():
        return {}

    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return {}
        code_field = next(
            (field for field in reader.fieldnames if field.lower().startswith("code")),
            None,
        )
        if code_field is None:
            raise ValueError(f"No population-code column found in {csv_path}")

        codes: dict[str, PopulationCode] = {}
        for row in reader:
            code = clean_cell(row.get(code_field))
            if not code:
                continue
            codes[code] = PopulationCode(
                code=code,
                species=clean_cell(row.get("Species")),
                population_name=clean_cell(row.get("Population Name")),
            )
    return codes


def infer_species_from_popcode(popcode: str) -> str:
    """Infer species from the main dataset popcode when the CSV is blank."""

    if popcode.startswith("CY_"):
        return "D. cymosa"
    if popcode.startswith("AB"):
        return "D. abramsii"
    if popcode:
        return "D. setchellii"
    return ""


def build_manifest(
    fastq_paths: Iterable[Path],
    population_codes: dict[str, PopulationCode],
) -> tuple[list[ManifestRow], list[ManifestIssue]]:
    """Build one manifest row per biological sample and collect issues."""

    records_by_sample: dict[tuple[str, str], list[FastqRecord]] = defaultdict(list)
    issues: list[ManifestIssue] = []

    for path in fastq_paths:
        try:
            record = parse_fastq_path(path)
        except ValueError as error:
            issues.append(
                ManifestIssue(
                    sample_id=path.stem,
                    batch=infer_batch(path),
                    issue_type="unparsed_fastq_name",
                    details=str(error),
                )
            )
            continue
        records_by_sample[(record.batch, record.sample_id)].append(record)

    rows: list[ManifestRow] = []
    for (batch, sample_id), records in sorted(records_by_sample.items()):
        first = records[0]
        r1_records = sorted(
            [record for record in records if record.read == "R1"],
            key=lambda record: record.path.as_posix(),
        )
        r2_records = sorted(
            [record for record in records if record.read == "R2"],
            key=lambda record: record.path.as_posix(),
        )

        pair_status = determine_pair_status(r1_records, r2_records)
        if pair_status != "complete":
            issues.append(
                ManifestIssue(
                    sample_id=sample_id,
                    batch=batch,
                    issue_type=pair_status,
                    details=f"R1 files: {len(r1_records)}; R2 files: {len(r2_records)}",
                )
            )

        population = population_codes.get(first.popcode)
        species = ""
        population_name = ""
        metadata_status = ""
        if first.naming_profile == "main_standard":
            species = population.species if population else ""
            if not species:
                species = infer_species_from_popcode(first.popcode)
            population_name = population.population_name if population else ""
            metadata_status = "resolved" if population else "popcode_not_in_csv"
        elif first.naming_profile.startswith("initial_"):
            metadata_status = "unresolved_initial_sample"
        else:
            metadata_status = "unrecognized_filename_profile"

        analysis_status, analysis_note = determine_analysis_status(pair_status)

        rows.append(
            ManifestRow(
                sample_id=sample_id,
                batch=batch,
                naming_profile=first.naming_profile,
                popcode=first.popcode,
                species=species,
                population_name=population_name,
                du_id=first.du_id,
                lp_id=first.lp_id,
                sequencing_samples=join_unique(
                    record.sequencing_sample for record in records
                ),
                lanes=join_unique(record.lane for record in records),
                r1_paths=join_paths(record.path for record in r1_records),
                r2_paths=join_paths(record.path for record in r2_records),
                r1_count=len(r1_records),
                r2_count=len(r2_records),
                pair_status=pair_status,
                metadata_status=metadata_status,
                analysis_status=analysis_status,
                analysis_note=analysis_note,
            )
        )

    return rows, sorted(issues, key=lambda issue: (issue.batch, issue.sample_id))


def determine_pair_status(
    r1_records: list[FastqRecord],
    r2_records: list[FastqRecord],
) -> str:
    if not r1_records and not r2_records:
        return "missing_R1_and_R2"
    if not r1_records:
        return "missing_R1"
    if not r2_records:
        return "missing_R2"
    if len(r1_records) != len(r2_records):
        return "uneven_read_counts"
    if len(r1_records) > 1:
        return "complete_multi_file"
    return "complete"


def determine_analysis_status(pair_status: str) -> tuple[str, str]:
    """Translate read-pair status into the primary alignment decision."""

    if pair_status == "complete":
        return (
            "include_primary_paired_end",
            "Use in the primary paired-end cpDNA/mtDNA alignment workflow.",
        )
    if pair_status in {"missing_R1", "missing_R2", "missing_R1_and_R2"}:
        return (
            "exclude_missing_mate",
            "Exclude from the primary paired-end cpDNA/mtDNA analysis because "
            "the mate FASTQ is absent after manual verification. If this sample "
            "is ever used for an individual single-end alignment, report it "
            "separately as a sensitivity check and do not mix it into the "
            "primary paired-end dataset.",
        )
    return (
        "review_before_primary_analysis",
        "Review before primary alignment because this sample does not have "
        "exactly one R1 and one R2 file.",
    )


def join_unique(values: Iterable[str]) -> str:
    return ";".join(sorted({value for value in values if value}))


def join_paths(paths: Iterable[Path]) -> str:
    return ";".join(path.as_posix() for path in paths)


def write_manifest_outputs(
    rows: list[ManifestRow],
    issues: list[ManifestIssue],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_dataclass_tsv(output_dir / "samples.tsv", rows, ManifestRow)
    write_dataclass_tsv(
        output_dir / "analysis_samples.tsv",
        primary_analysis_rows(rows),
        ManifestRow,
    )
    write_dataclass_tsv(
        output_dir / "excluded_samples.tsv",
        excluded_analysis_rows(rows),
        ManifestRow,
    )
    write_dataclass_tsv(output_dir / "pairing_report.tsv", issues, ManifestIssue)
    write_preflight_summary(output_dir / "preflight_summary.md", rows, issues)


def primary_analysis_rows(rows: list[ManifestRow]) -> list[ManifestRow]:
    return [
        row for row in rows if row.analysis_status == "include_primary_paired_end"
    ]


def excluded_analysis_rows(rows: list[ManifestRow]) -> list[ManifestRow]:
    return [
        row for row in rows if row.analysis_status != "include_primary_paired_end"
    ]


def write_dataclass_tsv(path: Path, rows: list[object], row_type: type[object]) -> None:
    fieldnames = list(row_type.__dataclass_fields__.keys())  # type: ignore[attr-defined]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: getattr(row, field) for field in fieldnames})


def write_preflight_summary(
    path: Path,
    rows: list[ManifestRow],
    issues: list[ManifestIssue],
) -> None:
    total = len(rows)
    complete = sum(1 for row in rows if row.pair_status == "complete")
    primary = sum(
        1 for row in rows if row.analysis_status == "include_primary_paired_end"
    )
    excluded = total - primary
    metadata_resolved = sum(1 for row in rows if row.metadata_status == "resolved")

    by_profile: dict[str, int] = defaultdict(int)
    by_batch: dict[str, int] = defaultdict(int)
    by_pair_status: dict[str, int] = defaultdict(int)
    by_metadata_status: dict[str, int] = defaultdict(int)
    by_analysis_status: dict[str, int] = defaultdict(int)

    for row in rows:
        by_profile[row.naming_profile] += 1
        by_batch[row.batch] += 1
        by_pair_status[row.pair_status] += 1
        by_metadata_status[row.metadata_status] += 1
        by_analysis_status[row.analysis_status] += 1

    lines = [
        "# Dudleya Organelle Alignment Preflight Summary",
        "",
        "This is step 1 of the cpDNA/mtDNA alignment pipeline. It validates",
        "sample naming, R1/R2 pairing, and population-code metadata before",
        "any read alignment is attempted.",
        "",
        "## Overall",
        "",
        f"- Samples discovered: {total}",
        f"- Samples with exactly one R1 and one R2: {complete}",
        f"- Samples in primary paired-end alignment set: {primary}",
        f"- Samples excluded from primary paired-end alignment: {excluded}",
        f"- Samples with resolved population metadata: {metadata_resolved}",
        f"- Issues reported: {len(issues)}",
        "",
        "## Samples By Sequencing Batch",
        "",
        *format_counts(by_batch),
        "",
        "## Samples By Naming Profile",
        "",
        *format_counts(by_profile),
        "",
        "## R1/R2 Pair Status",
        "",
        *format_counts(by_pair_status),
        "",
        "## Metadata Status",
        "",
        *format_counts(by_metadata_status),
        "",
        "## Primary Analysis Status",
        "",
        *format_counts(by_analysis_status),
        "",
        "## Missing-Mate Policy",
        "",
        "Samples without both mates are excluded from the primary paired-end",
        "cpDNA/mtDNA alignment. They remain documented in `samples.tsv`,",
        "`excluded_samples.tsv`, and `pairing_report.tsv`. If any missing-mate",
        "sample is ever aligned as an individual single-end case, that run must",
        "be reported separately as a sensitivity check and must not be mixed into",
        "the primary paired-end dataset.",
        "",
        "## Notes For The Next Pipeline Step",
        "",
        "- `analysis_samples.tsv` is the input table for primary paired-end",
        "  pilot read-to-reference alignment.",
        "- `excluded_samples.tsv` records samples excluded from the primary",
        "  alignment set and why.",
        "- `main_standard` samples can be used for population-level analyses when",
        "  their popcode appears in the population-code CSV.",
        "- `initial_du_dash` and `initial_du_lp` samples can be aligned, but should",
        "  remain metadata-unresolved until a manual lookup table is added.",
        "- No alignment, trimming, variant calling, or consensus generation happens",
        "  in this step.",
        "",
    ]

    path.write_text("\n".join(lines))


def format_counts(counts: dict[str, int]) -> list[str]:
    if not counts:
        return ["- none"]
    return [f"- {key}: {counts[key]}" for key in sorted(counts)]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build the Dudleya organelle FASTQ manifest and preflight reports."
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("genomicsDrive_data_dump"),
        help="Root folder containing downloaded QB3 FASTQ datasets.",
    )
    parser.add_argument(
        "--population-codes",
        type=Path,
        default=Path(
            "genomicsDrive_data_dump/QB3.Berkeley.251217/"
            "Dudleya DNAx - Population Codes.csv"
        ),
        help="CSV mapping popcodes to population/species names.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("dudleya_organelle_alignment_pipeline/results/00_manifest"),
        help="Folder for samples.tsv, pairing_report.tsv, and preflight_summary.md.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    fastq_paths = discover_fastq_paths(args.data_root)
    population_codes = load_population_codes(args.population_codes)
    rows, issues = build_manifest(fastq_paths, population_codes)
    write_manifest_outputs(rows, issues, args.output_dir)
    print(f"FASTQ files discovered: {len(fastq_paths)}")
    print(f"Samples written: {len(rows)}")
    print(f"Issues written: {len(issues)}")
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
