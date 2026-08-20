#!/usr/bin/env python3
"""Independently cross-check canonical pi and Hudson FST with scikit-allel."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import allel
import numpy as np
from organelle_pipeline.crosscheck import ratio_of_jointly_defined_components
from organelle_pipeline.paths import validate_run_id

BASE_CODES = np.array([ord("A"), ord("C"), ord("G"), ord("T")], dtype=np.uint8)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = ""
    parts: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                records[name] = "".join(parts)
            name = line[1:].split()[0]
            parts = []
        else:
            parts.append(line.strip().upper())
    if name:
        records[name] = "".join(parts)
    return records


def allele_counts(encoded: np.ndarray) -> np.ndarray:
    return np.stack([(encoded == code).sum(axis=0) for code in BASE_CODES], axis=1)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def equal_or_nan(left: float, right: float, tolerance: float = 1e-12) -> bool:
    return (math.isnan(left) and math.isnan(right)) or math.isclose(left, right, rel_tol=0, abs_tol=tolerance)


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    result_dir = root / "canonical_publication/results/popgen" / args.run_id
    provenance_dir = root / "canonical_publication/provenance/runs" / args.run_id
    output = result_dir / "trusted_scikit_allel_crosscheck.tsv"
    state_path = provenance_dir / "trusted_crosscheck.json"
    input_paths = [
        root / base / args.run_id / name
        for base, name in (
            ("canonical_publication/metadata/qc", "chloroplast_samples.tsv"),
            ("canonical_publication/metadata/qc", "mitochondria_samples.tsv"),
            ("canonical_publication/results/alignments", "chloroplast.callable_alignment.fa"),
            ("canonical_publication/results/alignments", "mitochondria.callable_alignment.fa"),
            ("canonical_publication/results/popgen", "chloroplast.population_summary.tsv"),
            ("canonical_publication/results/popgen", "mitochondria.population_summary.tsv"),
            ("canonical_publication/results/popgen", "chloroplast.pairwise_hudson_fst.tsv"),
            ("canonical_publication/results/popgen", "mitochondria.pairwise_hudson_fst.tsv"),
        )
    ]
    popgen_states = {organelle: provenance_dir / "popgen" / f"{organelle}.json" for organelle in ("chloroplast", "mitochondria")}
    input_paths.extend(popgen_states.values())
    upstream_fingerprints = {organelle: json.loads(path.read_text())["fingerprint"]["digest"] for organelle, path in popgen_states.items()}
    fingerprint_payload = {
        "implementation": "scikit-allel independent pi and Hudson FST",
        "command": "independent callable-site pi and Hudson ratio-of-sums FST cross-check",
        "git_commit": subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip(),
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "scikit_allel_version": allel.__version__,
        "python_executable_sha256": sha256_file(Path(sys.executable).resolve()),
        "script_sha256": sha256_file(Path(__file__)),
        "inputs": {path.relative_to(root).as_posix(): sha256_file(path) for path in input_paths},
        "upstream_fingerprints": upstream_fingerprints,
    }
    fingerprint = hashlib.sha256(json.dumps(fingerprint_payload, sort_keys=True).encode()).hexdigest()
    if args.resume and state_path.exists():
        saved = json.loads(state_path.read_text())
        if saved.get("fingerprint") != fingerprint:
            raise RuntimeError("Trusted cross-check state is stale")
        if not output.is_file() or sha256_file(output) != saved.get("output_sha256"):
            raise RuntimeError("Trusted cross-check output checksum mismatch")
        print("resume-valid trusted scikit-allel cross-check")
        return 0
    if state_path.exists() or output.exists():
        raise RuntimeError("Existing unvalidated trusted cross-check output")
    rows: list[dict[str, object]] = []
    failures: list[str] = []
    for organelle in ("chloroplast", "mitochondria"):
        metadata_rows = read_tsv(root / "canonical_publication/metadata/qc" / args.run_id / f"{organelle}_samples.tsv")
        population_by_sample = {row["sample_id"]: row["popcode"] for row in metadata_rows if row["popcode"]}
        records = read_fasta(root / "canonical_publication/results/alignments" / args.run_id / f"{organelle}.callable_alignment.fa")
        groups: dict[str, list[str]] = defaultdict(list)
        for sample, population in population_by_sample.items():
            groups[population].append(records[sample])
        canonical_pi = {
            row["population"]: float(row["nucleotide_diversity"]) for row in read_tsv(result_dir / f"{organelle}.population_summary.tsv")
        }
        for population, sequences in sorted(groups.items()):
            encoded = np.stack([np.frombuffer(sequence.encode(), dtype=np.uint8) for sequence in sequences])
            joint = np.isin(encoded, BASE_CODES).all(axis=0)
            if len(sequences) < 2 or not joint.any():
                trusted = math.nan
            else:
                trusted = float(np.mean(allel.mean_pairwise_difference(allele_counts(encoded[:, joint]))))
            passed = equal_or_nan(canonical_pi[population], trusted)
            rows.append(
                {
                    "organelle": organelle,
                    "statistic": "nucleotide_diversity",
                    "comparison": population,
                    "canonical": f"{canonical_pi[population]:.15g}",
                    "scikit_allel": f"{trusted:.15g}",
                    "absolute_difference": f"{abs(canonical_pi[population] - trusted):.15g}",
                    "match": "yes" if passed else "no",
                }
            )
            if not passed:
                failures.append(f"{organelle} pi {population}")
        for canonical in read_tsv(result_dir / f"{organelle}.pairwise_hudson_fst.tsv"):
            left, right = canonical["population_1"], canonical["population_2"]
            encoded_left = np.stack([np.frombuffer(sequence.encode(), dtype=np.uint8) for sequence in groups[left]])
            encoded_right = np.stack([np.frombuffer(sequence.encode(), dtype=np.uint8) for sequence in groups[right]])
            numerator, denominator = allel.hudson_fst(
                allele_counts(encoded_left),
                allele_counts(encoded_right),
            )
            trusted = ratio_of_jointly_defined_components(numerator, denominator)
            comparisons = (
                (
                    "hudson_numerator",
                    float(canonical["numerator"]),
                    trusted.numerator,
                    1e-9,
                ),
                (
                    "hudson_denominator",
                    float(canonical["denominator"]),
                    trusted.denominator,
                    1e-9,
                ),
                (
                    "hudson_callable_sites",
                    float(canonical["callable_sites_with_at_least_two_calls_per_population"]),
                    float(trusted.jointly_defined_sites),
                    0.0,
                ),
                ("hudson_fst", float(canonical["hudson_fst"]), trusted.ratio, 1e-12),
            )
            for statistic, canonical_value, trusted_value, tolerance in comparisons:
                passed = equal_or_nan(canonical_value, trusted_value, tolerance=tolerance)
                rows.append(
                    {
                        "organelle": organelle,
                        "statistic": statistic,
                        "comparison": f"{left}__{right}",
                        "canonical": f"{canonical_value:.15g}",
                        "scikit_allel": f"{trusted_value:.15g}",
                        "absolute_difference": f"{abs(canonical_value - trusted_value):.15g}",
                        "match": "yes" if passed else "no",
                    }
                )
                if not passed:
                    failures.append(f"{organelle} {statistic} {left}/{right}")
    with output.open("w", newline="") as handle:
        fields = ["organelle", "statistic", "comparison", "canonical", "scikit_allel", "absolute_difference", "match"]
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    state_path.write_text(
        json.dumps(
            {
                "status": "PASS" if not failures else "FAIL",
                "implementation": "scikit-allel",
                "version": allel.__version__,
                "fingerprint": fingerprint,
                "fingerprint_payload": fingerprint_payload,
                "comparison_count": len(rows),
                "failures": failures,
                "output": output.relative_to(root).as_posix(),
                "output_sha256": sha256_file(output),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    if failures:
        raise RuntimeError(f"Trusted population-genetics cross-check failed: {failures[:5]}")
    print(f"scikit-allel {allel.__version__}: {len(rows)} population-genetics comparisons PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
