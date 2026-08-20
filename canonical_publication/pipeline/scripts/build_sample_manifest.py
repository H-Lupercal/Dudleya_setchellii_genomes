#!/usr/bin/env python3
"""Regenerate canonical sample metadata from immutable raw filenames and CSV."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path

from organelle_pipeline.inventory import ACCEPTABLE_SOURCE_VALIDATION_STATUSES
from organelle_pipeline.metadata import (
    discover_samples,
    read_population_codes,
    write_sample_manifest,
)
from organelle_pipeline.paths import repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--population-codes", type=Path, required=True)
    parser.add_argument("--samples-output", type=Path, required=True)
    parser.add_argument("--populations-output", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    raw_root = root / repository_relative(args.raw_root, root)
    population_codes_path = root / repository_relative(args.population_codes, root)
    samples_output = root / repository_relative(args.samples_output, root)
    populations_output = root / repository_relative(args.populations_output, root)
    config_path = root / repository_relative(args.config, root)
    source_manifest = root / "canonical_publication/provenance/manifests/source_inputs.tsv"
    source_state = root / "canonical_publication/provenance/runs" / args.run_id / "source_validation.json"
    state_path = root / "canonical_publication/provenance/runs" / args.run_id / "metadata.json"
    state_path.parent.mkdir(parents=True, exist_ok=True)
    source_validation = json.loads(source_state.read_text())
    if source_validation.get("status") not in ACCEPTABLE_SOURCE_VALIDATION_STATUSES:
        raise RuntimeError("Sample metadata requires successful immutable-source validation")
    validate_saved_outputs(root, source_validation)
    source_fingerprint = source_validation["fingerprint"]["digest"]
    fingerprint = build_stage_fingerprint_from_hashes(
        "metadata",
        {
            **runtime_provenance(root, {"python": ("python", "--version")}),
            config_path.relative_to(root).as_posix(): sha256_file(config_path),
            source_manifest.relative_to(root).as_posix(): sha256_file(source_manifest),
            source_state.relative_to(root).as_posix(): sha256_file(source_state),
            population_codes_path.relative_to(root).as_posix(): sha256_file(population_codes_path),
        },
        {"immutable_sources": source_fingerprint},
        [
            "discover FASTQ pairs from immutable source paths; join population codes; "
            "materialize source-declared DUSE default; preserve metadata ambiguities; "
            "no manual exclusions"
        ],
    )
    if args.resume and state_path.exists():
        saved = json.loads(state_path.read_text())
        validate_resume(saved["fingerprint"]["digest"], fingerprint)
        for path, digest in saved["outputs"].items():
            if sha256_file(root / path) != digest:
                raise RuntimeError(f"Metadata output checksum mismatch: {path}")
        print("resume-valid metadata")
        return 0
    if state_path.exists():
        raise RuntimeError("Metadata state already exists; use --resume or a new run ID")
    populations = read_population_codes(population_codes_path)
    samples = discover_samples(raw_root, populations)
    write_sample_manifest(samples_output, samples, root)
    populations_output.parent.mkdir(parents=True, exist_ok=True)
    with populations_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["popcode", "species", "population_name"])
        for code in sorted(populations):
            record = populations[code]
            writer.writerow([record.code, record.species, record.population_name])
    ambiguity_output = root / "canonical_publication/metadata/qc" / args.run_id / "source_metadata_ambiguities.tsv"
    ambiguity_output.parent.mkdir(parents=True, exist_ok=True)
    with population_codes_path.open(newline="", encoding="utf-8-sig") as handle:
        source_rows = list(csv.DictReader(handle))
    source_fields = source_rows[0].keys() if source_rows else ()
    source_code_field = next(name for name in source_fields if name.strip().lower().startswith("code"))
    source_species_field = next(name for name in source_fields if name.strip().lower() == "species")
    blank_species_codes = {
        (row.get(source_code_field) or "").strip() for row in source_rows if not (row.get(source_species_field) or "").strip()
    }
    with ambiguity_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["popcode", "issue", "preserved_value"])
        for code, record in sorted(populations.items()):
            if " | " in record.population_name:
                writer.writerow([code, "conflicting source population labels", record.population_name])
            if code == "DUSE" and record.population_name.startswith("source-declared"):
                writer.writerow(
                    [
                        code,
                        "population code declared only in source column header",
                        record.population_name,
                    ]
                )
            if code in blank_species_codes and record.species:
                writer.writerow(
                    [
                        code,
                        "blank source species inferred from explicit CY_ code convention",
                        record.species,
                    ]
                )
    outputs = {
        path.resolve().relative_to(root).as_posix(): sha256_file(path) for path in (samples_output, populations_output, ambiguity_output)
    }
    state_path.write_text(
        json.dumps(
            {
                "status": "complete",
                "fingerprint": asdict(fingerprint),
                "outputs": outputs,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    complete = sum(sample.pair_status == "complete" for sample in samples)
    print(f"discovered {len(samples)} samples; {complete} complete pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
