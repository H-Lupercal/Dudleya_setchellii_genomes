"""Regenerate sample and population metadata from immutable source files."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

FASTQ_PATTERN = re.compile(
    r"^(?P<sample>.+)_S\d+_L\d+_R(?P<mate>[12])_\d+\.f(?:ast)?q\.gz$",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class PopulationCode:
    code: str
    species: str
    population_name: str


@dataclass(frozen=True)
class SampleRecord:
    sample_id: str
    popcode: str
    species: str
    population_name: str
    r1_paths: tuple[Path, ...]
    r2_paths: tuple[Path, ...]
    pair_status: str


def read_population_codes(path: Path | str) -> dict[str, PopulationCode]:
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Population code table has no header: {path}")
        code_column = next(
            (name for name in reader.fieldnames if name.strip().lower().startswith("code")),
            None,
        )
        species_column = next(
            (name for name in reader.fieldnames if name.strip().lower() == "species"),
            None,
        )
        population_column = next(
            (name for name in reader.fieldnames if name.strip().lower() == "population name"),
            None,
        )
        if code_column is None:
            raise ValueError(f"Population code table lacks a code column: {path}")
        result: dict[str, PopulationCode] = {}
        for row in reader:
            code = (row.get(code_column) or "").strip()
            if not code:
                continue
            species = (row.get(species_column) or "").strip() if species_column else ""
            if not species and code.startswith("CY_"):
                species = "D. cymosa"
            population_name = (row.get(population_column) or "").strip() if population_column else ""
            if code in result:
                existing = result[code]
                if existing.species and species and existing.species != species:
                    raise ValueError(f"Conflicting species labels for population code {code}")
                species = existing.species or species
                population_name = " | ".join(sorted({label for label in (existing.population_name, population_name) if label}))
            result[code] = PopulationCode(
                code=code,
                species=species,
                population_name=population_name,
            )
        if "= duse" in code_column.lower() and "DUSE" not in result:
            result["DUSE"] = PopulationCode(
                code="DUSE",
                species="D. setchellii",
                population_name="source-declared unprefixed DUSE group",
            )
    return result


def _population_for_sample(sample_id: str, population_codes: dict[str, PopulationCode]) -> PopulationCode | None:
    for code in sorted(population_codes, key=len, reverse=True):
        if sample_id == code or sample_id.startswith(f"{code}_"):
            return population_codes[code]
    # The source table's code-column header explicitly declares that numbered
    # DU IDs lacking a two-letter population prefix belong to DUSE.
    if "DUSE" in population_codes and re.match(r"^DU(?:-|\d)", sample_id):
        return population_codes["DUSE"]
    return None


def discover_samples(
    raw_root: Path | str,
    population_codes: dict[str, PopulationCode],
) -> list[SampleRecord]:
    """Discover every sample, retaining incomplete pairs for explicit reporting."""

    grouped: dict[str, dict[int, list[Path]]] = defaultdict(lambda: {1: [], 2: []})
    for path in sorted(Path(raw_root).rglob("*.f*q.gz")):
        match = FASTQ_PATTERN.match(path.name)
        if match is None:
            continue
        grouped[match.group("sample")][int(match.group("mate"))].append(path)
    records: list[SampleRecord] = []
    for sample_id in sorted(grouped):
        r1 = tuple(sorted(grouped[sample_id][1]))
        r2 = tuple(sorted(grouped[sample_id][2]))
        if r1 and r2 and len(r1) == len(r2):
            status = "complete"
        elif r1 and not r2:
            status = "missing_R2"
        elif r2 and not r1:
            status = "missing_R1"
        else:
            status = "unbalanced_lanes"
        population = _population_for_sample(sample_id, population_codes)
        records.append(
            SampleRecord(
                sample_id=sample_id,
                popcode=population.code if population else "",
                species=population.species if population else "",
                population_name=population.population_name if population else "",
                r1_paths=r1,
                r2_paths=r2,
                pair_status=status,
            )
        )
    return records


def write_sample_manifest(
    destination: Path | str,
    samples: list[SampleRecord],
    repository_root: Path | str,
) -> None:
    root = Path(repository_root).resolve()

    def relative(paths: tuple[Path, ...]) -> str:
        return ";".join(path.resolve().relative_to(root).as_posix() for path in paths)

    output = Path(destination)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            lineterminator="\n",
            fieldnames=[
                "sample_id",
                "popcode",
                "species",
                "population_name",
                "r1_paths",
                "r2_paths",
                "r1_count",
                "r2_count",
                "pair_status",
                "analysis_eligible",
            ],
        )
        writer.writeheader()
        for sample in samples:
            writer.writerow(
                {
                    "sample_id": sample.sample_id,
                    "popcode": sample.popcode,
                    "species": sample.species,
                    "population_name": sample.population_name,
                    "r1_paths": relative(sample.r1_paths),
                    "r2_paths": relative(sample.r2_paths),
                    "r1_count": len(sample.r1_paths),
                    "r2_count": len(sample.r2_paths),
                    "pair_status": sample.pair_status,
                    "analysis_eligible": "yes" if sample.pair_status == "complete" else "no",
                }
            )
