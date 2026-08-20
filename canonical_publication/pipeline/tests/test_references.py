from pathlib import Path

import pytest
from organelle_pipeline.references import (
    ReferenceValidationError,
    normalize_chloroplast,
    read_single_fasta,
    write_combined_reference,
)


def test_chloroplast_normalization_removes_configured_terminal_redundancy() -> None:
    normalized, removed = normalize_chloroplast("ACGTACGTAC", deduplicated_length=8)
    assert normalized == "ACGTACGT"
    assert removed == "AC"


def test_chloroplast_normalization_rejects_invalid_length() -> None:
    with pytest.raises(ReferenceValidationError):
        normalize_chloroplast("ACGT", deduplicated_length=5)


def test_combined_reference_has_stable_organelle_names(tmp_path: Path) -> None:
    output = tmp_path / "combined.fa"
    write_combined_reference(output, "ACGT", "TTAA")
    assert read_single_fasta(output, expected_records=2) == [
        ("chloroplast", "ACGT"),
        ("mitochondria", "TTAA"),
    ]
