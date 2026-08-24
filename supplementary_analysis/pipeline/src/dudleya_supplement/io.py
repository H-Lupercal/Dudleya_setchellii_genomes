"""Small deterministic tabular and JSON I/O helpers."""

from __future__ import annotations

import csv
import json
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

from .paths import assert_output_path


def read_tsv(path: Path | str) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path | str, rows: Iterable[Mapping[str, object]], fields: Sequence[str], repository_root: Path | str) -> None:
    output = assert_output_path(path, repository_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path | str, value: object, repository_root: Path | str) -> None:
    output = assert_output_path(path, repository_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
