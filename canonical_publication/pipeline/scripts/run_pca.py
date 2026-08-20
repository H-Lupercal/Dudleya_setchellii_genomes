#!/usr/bin/env python3
"""Run organelle-specific PCA on the canonical MAC>=2 haploid SNP sets."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from dataclasses import asdict
from pathlib import Path

import numpy as np
from organelle_pipeline.ordination import prepare_haploid_pca_matrix
from organelle_pipeline.paths import validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)
from sklearn.decomposition import PCA


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def read_genotypes(vcf: Path) -> tuple[list[str], np.ndarray]:
    samples = subprocess.run(
        ["bcftools", "query", "-l", str(vcf)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.splitlines()
    completed = subprocess.run(
        ["bcftools", "query", "-f", "%POS[\t%GT]\n", str(vcf)],
        capture_output=True,
        text=True,
        check=True,
    )
    markers = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")[1:]
        markers.append([float(value) if value in {"0", "1"} else np.nan for value in fields])
    if not markers:
        raise RuntimeError(f"No MAC>=2 markers available for PCA: {vcf}")
    return samples, np.asarray(markers, dtype=float).T


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    variant_dir = root / "canonical_publication/results/variants" / args.run_id
    metadata_dir = root / "canonical_publication/metadata/qc" / args.run_id
    output_dir = root / "canonical_publication/results/pca" / args.run_id
    state_dir = root / "canonical_publication/provenance/runs" / args.run_id / "pca"
    for directory in (output_dir, state_dir):
        directory.mkdir(parents=True, exist_ok=True)
    for organelle in ("chloroplast", "mitochondria"):
        vcf = variant_dir / f"{organelle}.mac2_ordination.vcf.gz"
        metadata_path = metadata_dir / f"{organelle}_samples.tsv"
        variant_state = json.loads(
            (root / "canonical_publication/provenance/runs" / args.run_id / "variants" / f"{organelle}.json").read_text()
        )
        fingerprint = build_stage_fingerprint_from_hashes(
            f"pca:{organelle}",
            {
                **runtime_provenance(
                    root,
                    {
                        "bcftools": ("bcftools", "--version"),
                        "numpy": ("python", "-c", "import numpy; print(numpy.__version__)"),
                        "scikit-learn": ("python", "-c", "import sklearn; print(sklearn.__version__)"),
                    },
                ),
                vcf.relative_to(root).as_posix(): sha256_file(vcf),
                metadata_path.relative_to(root).as_posix(): sha256_file(metadata_path),
            },
            {"variants": variant_state["fingerprint"]["digest"]},
            ["mean imputation; marker standardization; sklearn PCA"],
        )
        coordinates_path = output_dir / f"{organelle}.coordinates.tsv"
        variance_path = output_dir / f"{organelle}.variance.tsv"
        state_path = state_dir / f"{organelle}.json"
        if args.resume and state_path.exists():
            saved = json.loads(state_path.read_text())
            validate_resume(saved["fingerprint"]["digest"], fingerprint)
            validate_saved_outputs(root, saved)
            print(f"resume-valid pca {organelle}")
            continue
        if state_path.exists() or any(path.exists() for path in (coordinates_path, variance_path)):
            raise RuntimeError(f"Existing unvalidated PCA output for {organelle}")
        samples, genotypes = read_genotypes(vcf)
        matrix = prepare_haploid_pca_matrix(genotypes)
        component_count = min(10, matrix.shape[0] - 1, matrix.shape[1])
        if component_count < 2:
            raise RuntimeError(
                f"At least two PCA components require at least two MAC>=2 markers "
                f"and three samples for {organelle}; observed matrix {matrix.shape}"
            )
        model = PCA(n_components=component_count, svd_solver="full")
        coordinates = model.fit_transform(matrix)
        with metadata_path.open(newline="") as handle:
            metadata = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        if set(samples) != set(metadata):
            raise RuntimeError(f"PCA VCF/metadata sample mismatch for {organelle}")
        with coordinates_path.open("w", newline="") as handle:
            fields = ["sample_id", "popcode", *[f"PC{i}" for i in range(1, component_count + 1)]]
            writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fields)
            writer.writeheader()
            for sample, values in zip(samples, coordinates, strict=True):
                writer.writerow(
                    {
                        "sample_id": sample,
                        "popcode": metadata[sample]["popcode"],
                        **{f"PC{i + 1}": f"{value:.12g}" for i, value in enumerate(values)},
                    }
                )
        with variance_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["component", "explained_variance_ratio"])
            for index, value in enumerate(model.explained_variance_ratio_, 1):
                writer.writerow([f"PC{index}", f"{value:.12g}"])
        outputs = {path.relative_to(root).as_posix(): sha256_file(path) for path in (coordinates_path, variance_path)}
        state_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "organelle": organelle,
                    "sample_count": len(samples),
                    "marker_count": matrix.shape[1],
                    "fingerprint": asdict(fingerprint),
                    "outputs": outputs,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(f"completed {organelle} PCA")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
