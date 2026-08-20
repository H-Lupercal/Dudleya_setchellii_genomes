#!/usr/bin/env python3
"""Generate organelle-specific complete-case haplotypes and MST networks."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from dataclasses import asdict
from itertools import combinations
from pathlib import Path

import networkx as nx
from organelle_pipeline.haplotypes import summarize_haplotypes
from organelle_pipeline.paths import validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)
from organelle_pipeline.references import read_single_fasta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def distance(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right, strict=True))


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    alignment_dir = root / "canonical_publication/results/alignments" / args.run_id
    metadata_dir = root / "canonical_publication/metadata/qc" / args.run_id
    output_dir = root / "canonical_publication/results/haplotypes" / args.run_id
    state_dir = root / "canonical_publication/provenance/runs" / args.run_id / "haplotypes"
    for directory in (output_dir, state_dir):
        directory.mkdir(parents=True, exist_ok=True)
    for organelle in ("chloroplast", "mitochondria"):
        metadata_path = metadata_dir / f"{organelle}_samples.tsv"
        with metadata_path.open(newline="") as handle:
            metadata = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
        alignment = alignment_dir / f"{organelle}.callable_alignment.fa"
        records = dict(read_single_fasta(alignment, expected_records=len(metadata)))
        consensus_state = json.loads(
            (root / "canonical_publication/provenance/runs" / args.run_id / "consensus" / f"{organelle}.json").read_text()
        )
        fingerprint = build_stage_fingerprint_from_hashes(
            f"haplotypes:{organelle}",
            {
                **runtime_provenance(
                    root,
                    {
                        "networkx": ("python", "-c", "import networkx; print(networkx.__version__)"),
                        "python": ("python", "--version"),
                    },
                ),
                alignment.relative_to(root).as_posix(): sha256_file(alignment),
                metadata_path.relative_to(root).as_posix(): sha256_file(metadata_path),
            },
            {"consensus": consensus_state["fingerprint"]["digest"]},
            ["all accepted variable sites including singletons; ambiguous samples unassigned; haplotype minimum-spanning network"],
        )
        assignments_path = output_dir / f"{organelle}.sample_haplotypes.tsv"
        haplotypes_path = output_dir / f"{organelle}.haplotypes.tsv"
        positions_path = output_dir / f"{organelle}.haplotype_positions.tsv"
        edges_path = output_dir / f"{organelle}.network_edges.tsv"
        state_path = state_dir / f"{organelle}.json"
        if args.resume and state_path.exists():
            saved = json.loads(state_path.read_text())
            validate_resume(saved["fingerprint"]["digest"], fingerprint)
            validate_saved_outputs(root, saved)
            print(f"resume-valid haplotypes {organelle}")
            continue
        if state_path.exists() or any(path.exists() for path in (assignments_path, haplotypes_path, positions_path, edges_path)):
            raise RuntimeError(f"Existing unvalidated haplotype output for {organelle}")
        summary = summarize_haplotypes(records)
        with assignments_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample_id", "popcode", "haplotype"])
            for sample in sorted(records):
                writer.writerow([sample, metadata[sample]["popcode"], summary.sample_haplotypes[sample]])
        populations_by_haplotype = {
            haplotype: Counter(
                metadata[sample]["popcode"] or "unresolved"
                for sample, assigned in summary.sample_haplotypes.items()
                if assigned == haplotype
            )
            for haplotype in summary.counts
        }
        with haplotypes_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["haplotype", "sample_count", "population_counts", "sequence"])
            for haplotype in sorted(summary.counts):
                population_counts = ",".join(
                    f"{population}:{count}" for population, count in sorted(populations_by_haplotype[haplotype].items())
                )
                writer.writerow(
                    [
                        haplotype,
                        summary.counts[haplotype],
                        population_counts,
                        summary.sequences[haplotype],
                    ]
                )
        with positions_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["alignment_position_1based"])
            writer.writerows([[position + 1] for position in summary.positions])
        graph = nx.Graph()
        for haplotype in sorted(summary.counts):
            graph.add_node(haplotype)
        for left, right in combinations(sorted(summary.counts), 2):
            graph.add_edge(left, right, weight=distance(summary.sequences[left], summary.sequences[right]))
        tree = nx.minimum_spanning_tree(graph, weight="weight", algorithm="kruskal")
        with edges_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["haplotype_1", "haplotype_2", "mutational_distance"])
            for left, right, values in sorted(tree.edges(data=True)):
                writer.writerow([left, right, values["weight"]])
        outputs = {
            path.relative_to(root).as_posix(): sha256_file(path)
            for path in (
                assignments_path,
                haplotypes_path,
                positions_path,
                edges_path,
            )
        }
        state_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "organelle": organelle,
                    "haplotype_count": len(summary.counts),
                    "retained_high_confidence_variable_sites": len(summary.positions),
                    "ambiguous_sample_count": sum(label == "AMBIGUOUS" for label in summary.sample_haplotypes.values()),
                    "fingerprint": asdict(fingerprint),
                    "outputs": outputs,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(f"completed {organelle} haplotype network")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
