#!/usr/bin/env python3
"""Build supplementary sample-by-sample raw nucleotide-difference matrices."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from dataclasses import asdict
from itertools import combinations
from pathlib import Path

from organelle_pipeline.paths import validate_run_id
from organelle_pipeline.popgen import pack_sequence, packed_pairwise_distance
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)
from organelle_pipeline.references import read_single_fasta

ORGANELLES = ("chloroplast", "mitochondria")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def read_sample_ids(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sample_ids = [row["sample_id"] for row in rows]
    if len(sample_ids) != len(set(sample_ids)):
        raise RuntimeError(f"Duplicate sample ID in metadata: {path}")
    return sample_ids


def write_matrix(path: Path, sample_ids: list[str], values: list[list[int]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", *sample_ids])
        for sample_id, row in zip(sample_ids, values, strict=True):
            writer.writerow([sample_id, *row])


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    run_state = root / "canonical_publication/provenance/runs" / args.run_id
    alignment_dir = root / "canonical_publication/results/alignments" / args.run_id
    metadata_dir = root / "canonical_publication/metadata/qc" / args.run_id
    output_dir = root / "canonical_publication/results/supplement" / args.run_id / "pairwise_distances"
    state_dir = run_state / "distances"

    for organelle in ORGANELLES:
        alignment_path = alignment_dir / f"{organelle}.callable_alignment.fa"
        metadata_path = metadata_dir / f"{organelle}_samples.tsv"
        consensus_state_path = run_state / "consensus" / f"{organelle}.json"
        if not all(path.is_file() for path in (alignment_path, metadata_path, consensus_state_path)):
            raise RuntimeError(f"Missing canonical distance input for {organelle}")
        consensus_state = json.loads(consensus_state_path.read_text())
        validate_saved_outputs(root, consensus_state)
        declared_alignment = consensus_state.get("outputs", {}).get(alignment_path.relative_to(root).as_posix())
        if declared_alignment != sha256_file(alignment_path):
            raise RuntimeError(f"Consensus state does not declare the callable alignment for {organelle}")

        sample_ids = read_sample_ids(metadata_path)
        records = dict(read_single_fasta(alignment_path, expected_records=len(sample_ids)))
        if set(records) != set(sample_ids):
            raise RuntimeError(f"Alignment/metadata sample mismatch for {organelle}")
        declared = {
            **runtime_provenance(
                root,
                {
                    "numpy": ("python", "-c", "import numpy; print(numpy.__version__)"),
                    "python": ("python", "--version"),
                },
            ),
            alignment_path.relative_to(root).as_posix(): sha256_file(alignment_path),
            metadata_path.relative_to(root).as_posix(): sha256_file(metadata_path),
            consensus_state_path.relative_to(root).as_posix(): sha256_file(consensus_state_path),
        }
        fingerprint = build_stage_fingerprint_from_hashes(
            f"sample_distances:{organelle}",
            declared,
            {"consensus": consensus_state["fingerprint"]["digest"]},
            ["packed A/C/G/T bitsets; pair-specific callability; raw substitutions; p-distance"],
        )
        differences_path = output_dir / f"{organelle}.sample_pairwise_differences.tsv"
        callable_path = output_dir / f"{organelle}.sample_pairwise_callable_sites.tsv"
        long_path = output_dir / f"{organelle}.sample_pairwise_distances.tsv"
        outputs_to_write = (differences_path, callable_path, long_path)
        state_path = state_dir / f"{organelle}.json"
        if args.resume and state_path.exists():
            saved = json.loads(state_path.read_text())
            validate_resume(saved["fingerprint"]["digest"], fingerprint)
            validate_saved_outputs(root, saved)
            print(f"resume-valid sample distances {organelle}")
            continue
        if state_path.exists() or any(path.exists() for path in outputs_to_write):
            raise RuntimeError(f"Existing unvalidated sample-distance output for {organelle}")

        packed = [pack_sequence(records[sample_id]) for sample_id in sample_ids]
        sample_count = len(sample_ids)
        differences = [[0] * sample_count for _ in range(sample_count)]
        callable_sites = [[0] * sample_count for _ in range(sample_count)]
        for index, sequence in enumerate(packed):
            callable_sites[index][index] = sequence.callable_mask.bit_count()
        long_rows = []
        difference_values = []
        for left_index, right_index in combinations(range(sample_count), 2):
            result = packed_pairwise_distance(packed[left_index], packed[right_index])
            differences[left_index][right_index] = differences[right_index][left_index] = result.differences
            callable_sites[left_index][right_index] = callable_sites[right_index][left_index] = result.sites_compared
            difference_values.append(result.differences)
            long_rows.append(
                {
                    "organelle": organelle,
                    "sample_1": sample_ids[left_index],
                    "sample_2": sample_ids[right_index],
                    "differences": result.differences,
                    "sites_compared": result.sites_compared,
                    "p_distance": f"{result.p_distance:.12g}" if math.isfinite(result.p_distance) else "nan",
                }
            )

        output_dir.mkdir(parents=True, exist_ok=True)
        state_dir.mkdir(parents=True, exist_ok=True)
        write_matrix(differences_path, sample_ids, differences)
        write_matrix(callable_path, sample_ids, callable_sites)
        with long_path.open("w", newline="") as handle:
            fields = ["organelle", "sample_1", "sample_2", "differences", "sites_compared", "p_distance"]
            writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fields)
            writer.writeheader()
            writer.writerows(long_rows)
        outputs = {path.relative_to(root).as_posix(): sha256_file(path) for path in outputs_to_write}
        state_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "organelle": organelle,
                    "sample_count": sample_count,
                    "pairwise_comparison_count": len(long_rows),
                    "minimum_pairwise_differences": min(difference_values) if difference_values else 0,
                    "median_pairwise_differences": statistics.median(difference_values) if difference_values else 0,
                    "maximum_pairwise_differences": max(difference_values) if difference_values else 0,
                    "fingerprint": asdict(fingerprint),
                    "outputs": outputs,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(f"completed sample distances {organelle} for {sample_count} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
