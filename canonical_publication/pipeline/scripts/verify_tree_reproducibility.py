#!/usr/bin/env python3
"""Rerun fixed-seed trees and require reproducible strongly supported splits."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import tomllib
from dataclasses import asdict
from pathlib import Path

from organelle_pipeline.analysis import build_iqtree_command
from organelle_pipeline.paths import repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)
from organelle_pipeline.tree_reproducibility import compare_unrooted_trees


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    work = root / "canonical_publication/work" / args.run_id / "tree_reproducibility"
    log_dir = root / "canonical_publication/provenance/runs" / args.run_id / "logs/tree_reproducibility"
    state_path = root / "canonical_publication/provenance/runs" / args.run_id / "tree_reproducibility.json"
    report = root / "canonical_publication/results/trees" / args.run_id / "fixed_seed_reproducibility.tsv"
    tree_options = {
        "model": str(config["phylogeny"]["model"]),
        "sh_alrt_replicates": int(config["phylogeny"]["shalrt_replicates"]),
        "ultrafast_bootstrap_replicates": int(config["phylogeny"]["ultrafast_bootstrap_replicates"]),
        "bnni": bool(config["phylogeny"]["bnni"]),
    }
    alignment_dir = root / "canonical_publication/results/alignments" / args.run_id
    supplemental_dir = root / "canonical_publication/results/supplement" / args.run_id
    analyses = {
        "chloroplast": (
            alignment_dir / "chloroplast.callable_alignment.fa",
            None,
            root / "canonical_publication/results/trees" / args.run_id / "chloroplast.primary.treefile",
            int(config["phylogeny"]["cp_seed"]),
        ),
        "mitochondria": (
            alignment_dir / "mitochondria.callable_alignment.fa",
            None,
            root / "canonical_publication/results/trees" / args.run_id / "mitochondria.primary.treefile",
            int(config["phylogeny"]["mt_seed"]),
        ),
        "concatenated_supplementary": (
            supplemental_dir / "shared_partitioned_concatenated.fa",
            supplemental_dir / "shared_partitioned_concatenated.nex",
            supplemental_dir / "concatenated.partitioned.treefile",
            int(config["phylogeny"]["concatenated_seed"]),
        ),
    }
    declared = {
        **runtime_provenance(root, {"iqtree3": ("iqtree3", "--version")}),
        config_path.relative_to(root).as_posix(): sha256_file(config_path),
    }
    tree_state_paths = {
        "chloroplast": root / "canonical_publication/provenance/runs" / args.run_id / "trees/chloroplast.json",
        "mitochondria": root / "canonical_publication/provenance/runs" / args.run_id / "trees/mitochondria.json",
        "concatenated_supplementary": root / "canonical_publication/provenance/runs" / args.run_id / "trees/concatenated.json",
    }
    tree_states = {label: json.loads(path.read_text()) for label, path in tree_state_paths.items()}
    for path in tree_state_paths.values():
        declared[path.relative_to(root).as_posix()] = sha256_file(path)
    commands = []
    for label, (alignment, partition, canonical_tree, seed) in analyses.items():
        declared[alignment.relative_to(root).as_posix()] = sha256_file(alignment)
        declared[canonical_tree.relative_to(root).as_posix()] = sha256_file(canonical_tree)
        if partition is not None:
            declared[partition.relative_to(root).as_posix()] = sha256_file(partition)
        prefix = work / label / "replicate"
        commands.append(
            build_iqtree_command(
                alignment.relative_to(root),
                prefix.relative_to(root),
                seed=seed,
                partition_file=partition.relative_to(root) if partition else None,
                **tree_options,
            )
            + f" -nt {int(config['phylogeny']['threads'])}"
        )
    fingerprint = build_stage_fingerprint_from_hashes(
        "tree_reproducibility",
        declared,
        {label: state["fingerprint"]["digest"] for label, state in tree_states.items()},
        commands,
    )
    if args.resume and state_path.exists():
        saved = json.loads(state_path.read_text())
        validate_resume(saved["fingerprint"]["digest"], fingerprint)
        validate_saved_outputs(root, saved)
        print("resume-valid fixed-seed tree reproducibility")
        return 0
    if state_path.exists() or report.exists():
        raise RuntimeError("Existing unvalidated tree-reproducibility output")
    work.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    failures = []
    for command, (label, (_, _, canonical_tree, _)) in zip(commands, analyses.items(), strict=True):
        (work / label).mkdir(parents=True, exist_ok=True)
        replicate_tree = work / label / "replicate.treefile"
        replicate_internal_log = work / label / "replicate.log"
        if replicate_tree.exists():
            expected_log_command = f"Command: {command}"
            if not replicate_internal_log.exists() or expected_log_command not in replicate_internal_log.read_text():
                raise RuntimeError(f"Existing replicate command cannot be validated: {label}")
            print(f"validated existing fixed-seed replicate command for {label}")
        else:
            with (log_dir / f"{label}.log").open("w") as log:
                subprocess.run(
                    ["bash", "-o", "pipefail", "-c", command],
                    cwd=root,
                    env=os.environ.copy(),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
        canonical_hash = sha256_file(canonical_tree)
        replicate_hash = sha256_file(replicate_tree)
        byte_match = canonical_hash == replicate_hash
        comparison = compare_unrooted_trees(canonical_tree, replicate_tree)
        rows.append(
            {
                "analysis": label,
                "canonical_tree_sha256": canonical_hash,
                "replicate_tree_sha256": replicate_hash,
                "byte_identical": "yes" if byte_match else "no",
                "taxa_equal": "yes" if comparison.taxa_equal else "no",
                "canonical_internal_splits": comparison.canonical_internal_splits,
                "replicate_internal_splits": comparison.replicate_internal_splits,
                "full_unrooted_rf": comparison.full_unrooted_rf,
                "canonical_strong_splits": comparison.canonical_strong_splits,
                "replicate_strong_splits": comparison.replicate_strong_splits,
                "strong_split_symmetric_difference": comparison.strong_split_symmetric_difference,
                "max_shared_strong_branch_length_difference": f"{comparison.max_shared_strong_branch_length_difference:.12g}",
                "strong_topology_reproduced": "yes" if comparison.strong_topology_reproduced else "no",
            }
        )
        if not comparison.strong_topology_reproduced:
            failures.append(label)
    with report.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    state_path.write_text(
        json.dumps(
            {
                "status": "PASS" if not failures else "FAIL",
                "fingerprint": asdict(fingerprint),
                "failures": failures,
                "outputs": {report.relative_to(root).as_posix(): sha256_file(report)},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    if failures:
        raise RuntimeError(f"Fixed-seed tree reproducibility failed: {failures}")
    print("all fixed-seed strongly supported unrooted split sets reproduced")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
