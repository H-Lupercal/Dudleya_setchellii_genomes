"""Validation for supplementary sample-level nucleotide-distance outputs."""

from __future__ import annotations

import csv
import math
import statistics
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path


@dataclass(frozen=True)
class DistanceOutputSummary:
    sample_count: int
    pair_count: int
    minimum_differences: int
    median_differences: float
    maximum_differences: int


def _read_integer_matrix(path: Path, expected_samples: tuple[str, ...]) -> list[list[int]]:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] != ["sample_id", *expected_samples]:
        raise ValueError(f"Distance matrix header/sample order mismatch: {path}")
    if [row[0] for row in rows[1:]] != list(expected_samples):
        raise ValueError(f"Distance matrix row/sample order mismatch: {path}")
    if any(len(row) != len(expected_samples) + 1 for row in rows[1:]):
        raise ValueError(f"Distance matrix is not square: {path}")
    try:
        matrix = [[int(value) for value in row[1:]] for row in rows[1:]]
    except ValueError as error:
        raise ValueError(f"Distance matrix contains a non-integer value: {path}") from error
    if any(value < 0 for row in matrix for value in row):
        raise ValueError(f"Distance matrix contains a negative value: {path}")
    if any(matrix[left][right] != matrix[right][left] for left in range(len(matrix)) for right in range(len(matrix))):
        raise ValueError(f"Distance matrix must be symmetric: {path}")
    return matrix


def validate_pairwise_distance_outputs(
    organelle: str,
    expected_samples: tuple[str, ...],
    differences_path: Path,
    callable_sites_path: Path,
    long_form_path: Path,
    expected_callable_diagonal: tuple[int, ...] | None = None,
) -> DistanceOutputSummary:
    """Validate matrix/long-form agreement and return a raw-count summary."""

    differences = _read_integer_matrix(differences_path, expected_samples)
    callable_sites = _read_integer_matrix(callable_sites_path, expected_samples)
    sample_count = len(expected_samples)
    if any(differences[index][index] != 0 for index in range(sample_count)):
        raise ValueError(f"Difference matrix diagonal must be zero: {differences_path}")
    if expected_callable_diagonal is not None:
        if len(expected_callable_diagonal) != sample_count:
            raise ValueError("Expected callable diagonal/sample count mismatch")
        if any(callable_sites[index][index] != expected_callable_diagonal[index] for index in range(sample_count)):
            raise ValueError(f"Callable-site matrix diagonal mismatch: {callable_sites_path}")

    with long_form_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or ())
        rows = list(reader)
    required = {"organelle", "sample_1", "sample_2", "differences", "sites_compared", "p_distance"}
    if fields != required:
        raise ValueError(f"Invalid pairwise-distance columns: {long_form_path}")
    expected_pairs = list(combinations(range(sample_count), 2))
    if len(rows) != len(expected_pairs):
        raise ValueError(f"Pairwise-distance row count mismatch for {organelle}")
    values = []
    for row, (left, right) in zip(rows, expected_pairs, strict=True):
        expected_left, expected_right = expected_samples[left], expected_samples[right]
        if row["organelle"] != organelle or row["sample_1"] != expected_left or row["sample_2"] != expected_right:
            raise ValueError(f"Pairwise-distance sample order mismatch for {organelle}")
        observed_differences = int(row["differences"])
        observed_sites = int(row["sites_compared"])
        if observed_differences != differences[left][right] or observed_sites != callable_sites[left][right]:
            raise ValueError(f"Matrix/long-form distance mismatch for {expected_left}/{expected_right}")
        observed_p = float(row["p_distance"])
        expected_p = observed_differences / observed_sites if observed_sites else math.nan
        if not ((math.isnan(observed_p) and math.isnan(expected_p)) or math.isclose(observed_p, expected_p, rel_tol=1e-11, abs_tol=1e-15)):
            raise ValueError(f"Invalid p-distance for {expected_left}/{expected_right}")
        values.append(observed_differences)
    return DistanceOutputSummary(
        sample_count=sample_count,
        pair_count=len(values),
        minimum_differences=min(values) if values else 0,
        median_differences=statistics.median(values) if values else 0,
        maximum_differences=max(values) if values else 0,
    )
