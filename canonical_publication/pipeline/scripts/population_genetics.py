#!/usr/bin/env python3
"""Compute callable-site diversity, private variants, and signed Hudson FST."""

from __future__ import annotations

import argparse
import csv
import json
import math
import tomllib
from collections import defaultdict
from dataclasses import asdict
from itertools import combinations
from pathlib import Path

import numpy as np
from organelle_pipeline.paths import repository_relative, validate_run_id
from organelle_pipeline.popgen import (
    block_bootstrap_hudson_fst,
    callable_nucleotide_diversity,
    haplotype_diversity_from_assignments,
    hudson_fst,
    private_variant_sites_all,
)
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
)
from organelle_pipeline.references import read_single_fasta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numpy_pi(sequences: list[str]) -> tuple[int, int, int, float]:
    encoded = np.stack([np.frombuffer(sequence.encode(), dtype=np.uint8) for sequence in sequences])
    valid_codes = np.array([ord("A"), ord("C"), ord("G"), ord("T")], dtype=np.uint8)
    jointly_callable = np.isin(encoded, valid_codes).all(axis=0)
    differences = 0
    for left, right in combinations(encoded, 2):
        differences += int(((left != right) & jointly_callable).sum())
    joint_count = int(jointly_callable.sum())
    compared = joint_count * (len(sequences) * (len(sequences) - 1) // 2)
    return differences, compared, joint_count, differences / compared if compared else math.nan


def write_rows(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    block_size = int(config["population_genetics"]["bootstrap_block_size"])
    replicates = int(config["population_genetics"]["bootstrap_replicates"])
    seed = int(config["population_genetics"]["bootstrap_seed"])
    alignment_dir = root / "canonical_publication/results/alignments" / args.run_id
    metadata_dir = root / "canonical_publication/metadata/qc" / args.run_id
    result_dir = root / "canonical_publication/results/popgen" / args.run_id
    state_dir = root / "canonical_publication/provenance/runs" / args.run_id / "popgen"
    result_dir.mkdir(parents=True, exist_ok=True)
    state_dir.mkdir(parents=True, exist_ok=True)
    for organelle in ("chloroplast", "mitochondria"):
        metadata_path = metadata_dir / f"{organelle}_samples.tsv"
        metadata_rows = read_tsv(metadata_path)
        metadata = {row["sample_id"]: row for row in metadata_rows}
        alignment = alignment_dir / f"{organelle}.callable_alignment.fa"
        haplotype_assignments_path = root / "canonical_publication/results/haplotypes" / args.run_id / f"{organelle}.sample_haplotypes.tsv"
        reference_path = root / f"canonical_publication/references/selected/{organelle}.fa"
        reference = read_single_fasta(reference_path)[0][1]
        records = dict(read_single_fasta(alignment, expected_records=len(metadata)))
        if set(records) != set(metadata):
            raise RuntimeError(f"Alignment/metadata sample mismatch for {organelle}")
        consensus_state_path = root / "canonical_publication/provenance/runs" / args.run_id / "consensus" / f"{organelle}.json"
        consensus_state = json.loads(consensus_state_path.read_text())
        haplotype_state_path = root / "canonical_publication/provenance/runs" / args.run_id / "haplotypes" / f"{organelle}.json"
        haplotype_state = json.loads(haplotype_state_path.read_text())
        haplotype_rows = read_tsv(haplotype_assignments_path)
        haplotype_by_sample = {row["sample_id"]: row["haplotype"] for row in haplotype_rows}
        if len(haplotype_by_sample) != len(haplotype_rows) or set(haplotype_by_sample) != set(metadata):
            raise RuntimeError(f"Haplotype/metadata sample mismatch for {organelle}")
        if any(row["popcode"] != metadata[row["sample_id"]]["popcode"] for row in haplotype_rows):
            raise RuntimeError(f"Haplotype/metadata population mismatch for {organelle}")
        declared = {
            **runtime_provenance(
                root,
                {
                    "python": ("python", "--version"),
                    "numpy": ("python", "-c", "import numpy; print(numpy.__version__)"),
                },
            ),
            alignment.relative_to(root).as_posix(): sha256_file(alignment),
            haplotype_assignments_path.relative_to(root).as_posix(): sha256_file(haplotype_assignments_path),
            metadata_path.relative_to(root).as_posix(): sha256_file(metadata_path),
            reference_path.relative_to(root).as_posix(): sha256_file(reference_path),
            config_path.relative_to(root).as_posix(): sha256_file(config_path),
        }
        command = f"callable-site pi; Hudson ratio-of-sums FST; {block_size}-bp block bootstrap x{replicates}; seed={seed}"
        fingerprint = build_stage_fingerprint_from_hashes(
            f"popgen:{organelle}",
            declared,
            {
                "consensus": consensus_state["fingerprint"]["digest"],
                "haplotypes": haplotype_state["fingerprint"]["digest"],
            },
            [command],
        )
        population_output = result_dir / f"{organelle}.population_summary.tsv"
        pairwise_output = result_dir / f"{organelle}.pairwise_hudson_fst.tsv"
        crosscheck_output = result_dir / f"{organelle}.pi_crosscheck.tsv"
        state_path = state_dir / f"{organelle}.json"
        if args.resume and state_path.exists():
            saved = json.loads(state_path.read_text())
            validate_resume(saved["fingerprint"]["digest"], fingerprint)
            for path, digest in saved["outputs"].items():
                if sha256_file(root / path) != digest:
                    raise RuntimeError(f"Population-genetic output checksum mismatch: {path}")
            print(f"resume-valid popgen {organelle}")
            continue
        if state_path.exists() or any(path.exists() for path in (population_output, pairwise_output, crosscheck_output)):
            raise RuntimeError(f"Existing unvalidated population-genetics output for {organelle}")
        groups: dict[str, list[str]] = defaultdict(list)
        group_samples: dict[str, list[str]] = defaultdict(list)
        for sample, sequence in records.items():
            population = metadata[sample]["popcode"]
            if population:
                groups[population].append(sequence)
                group_samples[population].append(sample)
        observed_private = private_variant_sites_all(groups, reference)
        strict_private = private_variant_sites_all(
            groups,
            reference,
            require_joint_callability=True,
        )
        population_rows = []
        crosscheck_rows = []
        for population in sorted(groups):
            sequences = groups[population]
            pi = callable_nucleotide_diversity(sequences)
            independent = numpy_pi(sequences)
            if (
                pi.differences != independent[0]
                or pi.compared_sites != independent[1]
                or pi.jointly_callable_sites != independent[2]
                or not ((math.isnan(pi.pi) and math.isnan(independent[3])) or math.isclose(pi.pi, independent[3], rel_tol=0, abs_tol=1e-15))
            ):
                raise RuntimeError(f"Independent pi cross-check failed for {population}")
            haplotype = haplotype_diversity_from_assignments([haplotype_by_sample[sample] for sample in group_samples[population]])
            population_rows.append(
                {
                    "organelle": organelle,
                    "population": population,
                    "sample_count": len(sequences),
                    "sample_ids": ",".join(group_samples[population]),
                    "pairwise_differences": pi.differences,
                    "jointly_callable_sites": pi.jointly_callable_sites,
                    "pairwise_callable_sites": pi.compared_sites,
                    "nucleotide_diversity": f"{pi.pi:.12g}",
                    "organelle_global_haplotype_sites": haplotype_state["retained_high_confidence_variable_sites"],
                    "haplotype_count": haplotype.haplotype_count,
                    "haplotype_diversity": f"{haplotype.diversity:.12g}",
                    "haplotype_assigned_samples": haplotype.assigned_samples,
                    "haplotype_ambiguous_samples": haplotype.ambiguous_samples,
                    "observed_private_nonreference_sites_including_singletons": len(observed_private[population]),
                    "strict_jointly_callable_private_nonreference_sites": len(strict_private[population]),
                }
            )
            crosscheck_rows.append(
                {
                    "population": population,
                    "implementation": "numpy_pairwise_independent",
                    "differences": independent[0],
                    "compared_sites": independent[1],
                    "jointly_callable_sites": independent[2],
                    "pi": f"{independent[3]:.12g}",
                    "exact_match": "yes",
                }
            )
        pairwise_rows = []
        for pair_index, (left, right) in enumerate(combinations(sorted(groups), 2)):
            estimate = hudson_fst(groups[left], groups[right])
            lower, upper = block_bootstrap_hudson_fst(
                groups[left],
                groups[right],
                block_size=block_size,
                replicates=replicates,
                seed=seed + pair_index,
            )
            pairwise_rows.append(
                {
                    "organelle": organelle,
                    "population_1": left,
                    "population_2": right,
                    "n_population_1": len(groups[left]),
                    "n_population_2": len(groups[right]),
                    "numerator": f"{estimate.numerator:.12g}",
                    "denominator": f"{estimate.denominator:.12g}",
                    "callable_sites_with_at_least_two_calls_per_population": estimate.callable_sites,
                    "hudson_fst": f"{estimate.fst:.12g}",
                    "bootstrap_ci_2.5": f"{lower:.12g}",
                    "bootstrap_ci_97.5": f"{upper:.12g}",
                    "bootstrap_block_bp": block_size,
                    "bootstrap_replicates": replicates,
                    "bootstrap_seed": seed + pair_index,
                }
            )
        expected_pairs = len(groups) * (len(groups) - 1) // 2
        if len(pairwise_rows) != expected_pairs:
            raise RuntimeError("Population pair count invariant failed")
        write_rows(population_output, population_rows, list(population_rows[0]))
        write_rows(pairwise_output, pairwise_rows, list(pairwise_rows[0]))
        write_rows(crosscheck_output, crosscheck_rows, list(crosscheck_rows[0]))
        outputs = {path.relative_to(root).as_posix(): sha256_file(path) for path in (population_output, pairwise_output, crosscheck_output)}
        state_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "organelle": organelle,
                    "population_count": len(groups),
                    "pairwise_comparison_count": len(pairwise_rows),
                    "fingerprint": asdict(fingerprint),
                    "outputs": outputs,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(f"computed {organelle} population genetics for {len(groups)} populations")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
