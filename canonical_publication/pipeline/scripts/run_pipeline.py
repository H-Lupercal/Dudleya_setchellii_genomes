#!/usr/bin/env python3
"""Canonical dependency-ordered publication pipeline entrypoint."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tomllib
from pathlib import Path

from organelle_pipeline.configuration import validate_publication_config
from organelle_pipeline.paths import repository_relative, validate_run_id

STAGES = (
    "source_validation",
    "references",
    "metadata",
    "mapping",
    "mapping_provenance",
    "qc",
    "variants",
    "consensus",
    "pca",
    "haplotypes",
    "popgen",
    "crosscheck",
    "trees",
    "treecheck",
    "admixture",
    "figures",
    "reports",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--from-stage", choices=STAGES, default=STAGES[0])
    parser.add_argument("--until-stage", choices=STAGES, default=STAGES[-1])
    parser.add_argument("--validation-python")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = Path(__file__).resolve().parents[3]
    scripts = root / "canonical_publication/pipeline/scripts"
    validation_python = (
        args.validation_python or os.environ.get("SCIKIT_ALLEL_PYTHON") or str(root / ".tools/scikit-allel-validation/bin/python")
    )
    config = root / repository_relative(args.config, root)
    if not config.exists():
        raise SystemExit(f"Configuration not found: {config}")
    publication_config = tomllib.loads(config.read_text())
    validate_publication_config(publication_config)
    execution = publication_config["execution"]
    population_codes = root / "source_data/raw_reads/genomicsDrive_data_dump/QB3.Berkeley.251217" / "Dudleya DNAx - Population Codes.csv"
    resume = ["--resume"] if args.resume else []
    python = sys.executable
    commands = {
        "source_validation": [
            python,
            str(scripts / "verify_source_md5.py"),
            "--repository-root",
            str(root),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "references": [
            python,
            str(scripts / "prepare_references.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "metadata": [
            python,
            str(scripts / "build_sample_manifest.py"),
            "--repository-root",
            str(root),
            "--raw-root",
            str(root / "source_data/raw_reads/genomicsDrive_data_dump"),
            "--population-codes",
            str(population_codes),
            "--samples-output",
            str(root / "canonical_publication/metadata/samples/samples.tsv"),
            "--populations-output",
            str(root / "canonical_publication/metadata/populations/populations.tsv"),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "mapping": [
            python,
            str(scripts / "map_samples.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            "--jobs",
            str(execution["mapping_jobs"]),
            "--threads-per-job",
            str(execution["mapping_threads_per_job"]),
            *resume,
        ],
        "mapping_provenance": [
            python,
            str(scripts / "finalize_mapping_provenance.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            "--threads-per-job",
            str(execution["mapping_threads_per_job"]),
            *resume,
        ],
        "qc": [
            python,
            str(scripts / "compute_qc.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "variants": [
            python,
            str(scripts / "call_variants.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "consensus": [
            python,
            str(scripts / "build_consensus.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "pca": [
            python,
            str(scripts / "run_pca.py"),
            "--repository-root",
            str(root),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "haplotypes": [
            python,
            str(scripts / "run_haplotypes.py"),
            "--repository-root",
            str(root),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "popgen": [
            python,
            str(scripts / "population_genetics.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "crosscheck": [
            validation_python,
            str(scripts / "crosscheck_scikit_allel.py"),
            "--repository-root",
            str(root),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "trees": [
            python,
            str(scripts / "run_trees.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "treecheck": [
            python,
            str(scripts / "verify_tree_reproducibility.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "admixture": [
            python,
            str(scripts / "run_admixture.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            "--jobs",
            str(execution["admixture_jobs"]),
            "--threads-per-job",
            str(execution["admixture_threads_per_job"]),
            *resume,
        ],
        "figures": [
            python,
            str(scripts / "render_figures.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            *resume,
        ],
        "reports": [
            python,
            str(scripts / "build_reports.py"),
            "--repository-root",
            str(root),
            "--config",
            str(config),
            "--run-id",
            args.run_id,
            *resume,
        ],
    }
    start = STAGES.index(args.from_stage)
    end = STAGES.index(args.until_stage)
    if start > end:
        raise SystemExit("--from-stage must not follow --until-stage")
    if start and not args.resume:
        raise SystemExit("--from-stage requires --resume because all preceding stages must be fingerprint-validated before downstream work")
    execution_start = 0 if start else start
    environment = os.environ.copy()
    source_path = str(root / "canonical_publication/pipeline/src")
    environment["PYTHONPATH"] = source_path + (os.pathsep + environment["PYTHONPATH"] if environment.get("PYTHONPATH") else "")
    tool_path = str(root / ".tools/bioconda-env/bin")
    environment["PATH"] = tool_path + os.pathsep + environment["PATH"]
    for stage in STAGES[execution_start : end + 1]:
        command = commands[stage]
        role = "UPSTREAM VALIDATION" if STAGES.index(stage) < start else "STAGE"
        print(f"{role} {stage}: {' '.join(command)}", flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=root, env=environment, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
