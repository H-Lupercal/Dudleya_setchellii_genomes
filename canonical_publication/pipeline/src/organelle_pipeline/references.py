"""Reference FASTA validation and deterministic normalization."""

from __future__ import annotations

from pathlib import Path

VALID_BASES = frozenset("ACGTN")


class ReferenceValidationError(ValueError):
    """Raised when a candidate reference violates canonical assumptions."""


def read_single_fasta(path: Path | str, expected_records: int = 1) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    for raw_line in Path(path).read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts).upper()))
            header = line[1:].split()[0]
            parts = []
        elif header is None:
            raise ReferenceValidationError(f"Sequence appears before FASTA header: {path}")
        else:
            parts.append(line)
    if header is not None:
        records.append((header, "".join(parts).upper()))
    if len(records) != expected_records:
        raise ReferenceValidationError(f"Expected {expected_records} FASTA records in {path}, found {len(records)}")
    for record_header, sequence in records:
        invalid = set(sequence) - VALID_BASES
        if not sequence or invalid:
            raise ReferenceValidationError(f"Invalid sequence for {record_header}: {sorted(invalid)}")
    return records


def normalize_chloroplast(
    raw_sequence: str,
    deduplicated_length: int,
) -> tuple[str, str]:
    """Split a candidate at the externally validated circularization boundary."""

    sequence = raw_sequence.upper()
    if deduplicated_length <= 0 or deduplicated_length > len(sequence):
        raise ReferenceValidationError(f"Invalid deduplicated length {deduplicated_length} for {len(sequence)} bp")
    invalid = set(sequence) - VALID_BASES
    if invalid:
        raise ReferenceValidationError(f"Invalid chloroplast bases: {sorted(invalid)}")
    return sequence[:deduplicated_length], sequence[deduplicated_length:]


def write_fasta(path: Path | str, records: list[tuple[str, str]]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def write_combined_reference(path: Path | str, chloroplast: str, mitochondria: str) -> None:
    write_fasta(path, [("chloroplast", chloroplast), ("mitochondria", mitochondria)])
