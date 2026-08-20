#!/usr/bin/env python3
"""Infer primary organelle trees and a supplementary partitioned concatenation."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import tomllib
from dataclasses import asdict
from pathlib import Path

from Bio import Phylo
from organelle_pipeline.analysis import (
    IQTreeSupport,
    alignment_callability_counts,
    alignment_site_counts,
    build_iqtree_command,
    is_strong_iqtree_support,
    parse_iqtree_support,
)
from organelle_pipeline.paths import repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)
from organelle_pipeline.references import read_single_fasta, write_fasta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def read_sample_ids(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        return [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")]


def validate_tree_taxa(tree_path: Path, expected_taxa: set[str]) -> None:
    tree = Phylo.read(tree_path, "newick")
    observed = [terminal.name for terminal in tree.get_terminals()]
    if len(observed) != len(set(observed)):
        raise RuntimeError(f"Tree contains duplicate terminal IDs: {tree_path}")
    if set(observed) != expected_taxa:
        missing = sorted(expected_taxa - set(observed))
        extra = sorted(set(observed) - expected_taxa)
        raise RuntimeError(f"Tree/QC sample mismatch for {tree_path}: missing={missing[:5]}, extra={extra[:5]}")


def canonical_splits(tree_path: Path, taxa: set[str]) -> dict[frozenset[str], IQTreeSupport]:
    tree = Phylo.read(tree_path, "newick")
    terminal_names = [terminal.name for terminal in tree.get_terminals()]
    if len(set(terminal_names)) != len(terminal_names):
        raise ValueError(f"Tree contains duplicate terminal names: {tree_path}")
    if not taxa <= set(terminal_names):
        raise ValueError(f"Tree is missing requested shared taxa: {tree_path}")
    splits = {}
    for clade in tree.get_nonterminals(order="postorder"):
        side = {terminal.name for terminal in clade.get_terminals()} & taxa
        if len(side) < 2 or len(taxa - side) < 2:
            continue
        complement = taxa - side
        canonical = side if (len(side), sorted(side)) <= (len(complement), sorted(complement)) else complement
        try:
            support = parse_iqtree_support(clade.name, clade.confidence)
        except ValueError:
            continue
        key = frozenset(canonical)
        previous = splits.get(key)
        if previous is None or (
            min(support.sh_alrt, support.ultrafast_bootstrap),
            support.sh_alrt,
            support.ultrafast_bootstrap,
        ) > (
            min(previous.sh_alrt, previous.ultrafast_bootstrap),
            previous.sh_alrt,
            previous.ultrafast_bootstrap,
        ):
            # Multiple full-tree edges can induce the same split after
            # restricting to shared taxa. Retain the strongest actual edge;
            # never transfer a label onto a synthetic pruned node.
            splits[key] = support
    return splits


def incompatible(left: frozenset[str], right: frozenset[str], taxa: set[str]) -> bool:
    return all(
        (
            left & right,
            left & (taxa - right),
            (taxa - left) & right,
            (taxa - left) & (taxa - right),
        )
    )


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    alignment_dir = root / "canonical_publication/results/alignments" / args.run_id
    metadata_dir = root / "canonical_publication/metadata/qc" / args.run_id
    output_dir = root / "canonical_publication/results/trees" / args.run_id
    supplemental_dir = root / "canonical_publication/results/supplement" / args.run_id
    state_dir = root / "canonical_publication/provenance/runs" / args.run_id / "trees"
    log_dir = root / "canonical_publication/provenance/runs" / args.run_id / "logs/trees"
    for directory in (output_dir, supplemental_dir, state_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)
    tree_options = {
        "model": str(config["phylogeny"]["model"]),
        "sh_alrt_replicates": int(config["phylogeny"]["shalrt_replicates"]),
        "ultrafast_bootstrap_replicates": int(config["phylogeny"]["ultrafast_bootstrap_replicates"]),
        "bnni": bool(config["phylogeny"]["bnni"]),
    }
    minimum_strong_sh_alrt = float(config["phylogeny"]["strong_conflict_minimum_sh_alrt"])
    minimum_strong_ufboot = float(config["phylogeny"]["strong_conflict_minimum_ultrafast_bootstrap"])
    treefiles = {}
    for organelle, seed_key in (("chloroplast", "cp_seed"), ("mitochondria", "mt_seed")):
        alignment = alignment_dir / f"{organelle}.callable_alignment.fa"
        expected_taxa = set(read_sample_ids(metadata_dir / f"{organelle}_samples.tsv"))
        consensus_state = json.loads(
            (root / "canonical_publication/provenance/runs" / args.run_id / "consensus" / f"{organelle}.json").read_text()
        )
        prefix = output_dir / f"{organelle}.primary"
        command = (
            build_iqtree_command(
                alignment.relative_to(root),
                prefix.relative_to(root),
                seed=int(config["phylogeny"][seed_key]),
                **tree_options,
            )
            + f" -nt {int(config['phylogeny']['threads'])}"
        )
        fingerprint = build_stage_fingerprint_from_hashes(
            f"tree:{organelle}",
            {
                **runtime_provenance(
                    root,
                    {
                        "biopython": ("python", "-c", "import Bio; print(Bio.__version__)"),
                        "iqtree3": ("iqtree3", "--version"),
                    },
                ),
                config_path.relative_to(root).as_posix(): sha256_file(config_path),
                alignment.relative_to(root).as_posix(): sha256_file(alignment),
            },
            {"consensus": consensus_state["fingerprint"]["digest"]},
            [command],
        )
        state_path = state_dir / f"{organelle}.json"
        treefile = Path(f"{prefix}.treefile")
        if args.resume and state_path.exists():
            saved = json.loads(state_path.read_text())
            validate_resume(saved["fingerprint"]["digest"], fingerprint)
            validate_saved_outputs(root, saved)
            validate_tree_taxa(treefile, expected_taxa)
            treefiles[organelle] = treefile
            print(f"resume-valid tree {organelle}")
            continue
        if state_path.exists() or treefile.exists():
            raise RuntimeError(f"Existing unvalidated tree output for {organelle}")
        log_path = log_dir / f"{organelle}.log"
        with log_path.open("w") as log:
            subprocess.run(
                ["bash", "-o", "pipefail", "-c", command],
                cwd=root,
                env=os.environ.copy(),
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        outputs = {
            path.relative_to(root).as_posix(): sha256_file(path) for path in (treefile, Path(f"{prefix}.contree"), Path(f"{prefix}.iqtree"))
        }
        validate_tree_taxa(treefile, expected_taxa)
        state_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "organelle": organelle,
                    "sample_count": len(expected_taxa),
                    "rooting": "unrooted",
                    "fingerprint": asdict(fingerprint),
                    "outputs": outputs,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        treefiles[organelle] = treefile
        print(f"completed primary {organelle} tree")

    shared = read_sample_ids(metadata_dir / "shared_samples.tsv")
    if len(shared) < 4:
        raise RuntimeError("Supplementary unrooted concatenated tree requires at least four shared samples")
    shared_set = set(shared)
    cp_records = dict(
        read_single_fasta(
            alignment_dir / "chloroplast.callable_alignment.fa",
            expected_records=len(read_sample_ids(metadata_dir / "chloroplast_samples.tsv")),
        )
    )
    mt_records = dict(
        read_single_fasta(
            alignment_dir / "mitochondria.callable_alignment.fa",
            expected_records=len(read_sample_ids(metadata_dir / "mitochondria_samples.tsv")),
        )
    )
    if not shared_set <= set(cp_records) or not shared_set <= set(mt_records):
        raise RuntimeError("Shared sample set is not the organelle intersection")
    concatenated = supplemental_dir / "shared_partitioned_concatenated.fa"
    cp_length = len(next(iter(cp_records.values())))
    mt_length = len(next(iter(mt_records.values())))
    partition = supplemental_dir / "shared_partitioned_concatenated.nex"
    contribution = supplemental_dir / "concatenated_site_contribution.tsv"
    cp_shared = {sample: cp_records[sample] for sample in shared}
    mt_shared = {sample: mt_records[sample] for sample in shared}
    concat_state_path = state_dir / "concatenated.json"
    if not (args.resume and concat_state_path.exists()):
        if concat_state_path.exists() or any(path.exists() for path in (concatenated, partition, contribution)):
            raise RuntimeError("Existing unvalidated supplementary concatenation input")
        write_fasta(
            concatenated,
            [(sample, cp_records[sample] + mt_records[sample]) for sample in shared],
        )
        partition.write_text(
            "#nexus\nbegin sets;\n"
            f"  charset chloroplast = 1-{cp_length};\n"
            f"  charset mitochondria = {cp_length + 1}-{cp_length + mt_length};\n"
            "end;\n"
        )
        cp_variable, cp_parsimony = alignment_site_counts(cp_shared)
        mt_variable, mt_parsimony = alignment_site_counts(mt_shared)
        cp_callability = alignment_callability_counts(cp_shared)
        mt_callability = alignment_callability_counts(mt_shared)
        total_length = cp_callability.coordinate_span_sites + mt_callability.coordinate_span_sites
        total_two_calls = cp_callability.sites_with_at_least_two_callable_samples + mt_callability.sites_with_at_least_two_callable_samples
        total_joint = cp_callability.jointly_callable_sites + mt_callability.jointly_callable_sites
        total_variable = cp_variable + mt_variable
        total_parsimony = cp_parsimony + mt_parsimony
        with contribution.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "partition",
                    "coordinate_span_sites",
                    "coordinate_span_fraction",
                    "sites_with_at_least_two_callable_shared_samples",
                    "at_least_two_callable_site_fraction",
                    "jointly_callable_shared_sample_sites",
                    "jointly_callable_site_fraction",
                    "variable_sites_including_singletons",
                    "variable_site_fraction",
                    "parsimony_informative_sites",
                    "parsimony_informative_site_fraction",
                ]
            )
            writer.writerow(
                [
                    "chloroplast",
                    cp_callability.coordinate_span_sites,
                    f"{cp_callability.coordinate_span_sites / total_length:.12g}",
                    cp_callability.sites_with_at_least_two_callable_samples,
                    f"{cp_callability.sites_with_at_least_two_callable_samples / total_two_calls:.12g}" if total_two_calls else "nan",
                    cp_callability.jointly_callable_sites,
                    f"{cp_callability.jointly_callable_sites / total_joint:.12g}" if total_joint else "nan",
                    cp_variable,
                    f"{cp_variable / total_variable:.12g}" if total_variable else "nan",
                    cp_parsimony,
                    f"{cp_parsimony / total_parsimony:.12g}" if total_parsimony else "nan",
                ]
            )
            writer.writerow(
                [
                    "mitochondria",
                    mt_callability.coordinate_span_sites,
                    f"{mt_callability.coordinate_span_sites / total_length:.12g}",
                    mt_callability.sites_with_at_least_two_callable_samples,
                    f"{mt_callability.sites_with_at_least_two_callable_samples / total_two_calls:.12g}" if total_two_calls else "nan",
                    mt_callability.jointly_callable_sites,
                    f"{mt_callability.jointly_callable_sites / total_joint:.12g}" if total_joint else "nan",
                    mt_variable,
                    f"{mt_variable / total_variable:.12g}" if total_variable else "nan",
                    mt_parsimony,
                    f"{mt_parsimony / total_parsimony:.12g}" if total_parsimony else "nan",
                ]
            )
    concat_prefix = supplemental_dir / "concatenated.partitioned"
    concat_command = (
        build_iqtree_command(
            concatenated.relative_to(root),
            concat_prefix.relative_to(root),
            seed=int(config["phylogeny"]["concatenated_seed"]),
            partition_file=partition.relative_to(root),
            **tree_options,
        )
        + f" -nt {int(config['phylogeny']['threads'])}"
    )
    concat_fingerprint = build_stage_fingerprint_from_hashes(
        "tree:concatenated_supplementary",
        {
            **runtime_provenance(
                root,
                {
                    "biopython": ("python", "-c", "import Bio; print(Bio.__version__)"),
                    "iqtree3": ("iqtree3", "--version"),
                },
            ),
            config_path.relative_to(root).as_posix(): sha256_file(config_path),
            concatenated.relative_to(root).as_posix(): sha256_file(concatenated),
            partition.relative_to(root).as_posix(): sha256_file(partition),
        },
        {
            "chloroplast": json.loads((state_dir / "chloroplast.json").read_text())["fingerprint"]["digest"],
            "mitochondria": json.loads((state_dir / "mitochondria.json").read_text())["fingerprint"]["digest"],
        },
        [concat_command],
    )
    concat_tree = Path(f"{concat_prefix}.treefile")
    if not (args.resume and concat_state_path.exists()):
        if concat_state_path.exists() or concat_tree.exists():
            raise RuntimeError("Existing unvalidated supplementary concatenated tree")
        with (log_dir / "concatenated.log").open("w") as log:
            subprocess.run(
                ["bash", "-o", "pipefail", "-c", concat_command],
                cwd=root,
                env=os.environ.copy(),
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        concat_state_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "role": "supplementary",
                    "sample_count": len(shared),
                    "fingerprint": asdict(concat_fingerprint),
                    "outputs": {
                        path.relative_to(root).as_posix(): sha256_file(path)
                        for path in (
                            concat_tree,
                            Path(f"{concat_prefix}.contree"),
                            Path(f"{concat_prefix}.iqtree"),
                            concatenated,
                            partition,
                            contribution,
                        )
                    },
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
    else:
        saved_concat = json.loads(concat_state_path.read_text())
        validate_resume(saved_concat["fingerprint"]["digest"], concat_fingerprint)
        validate_saved_outputs(root, saved_concat)
    validate_tree_taxa(concat_tree, shared_set)

    cp_splits = canonical_splits(treefiles["chloroplast"], shared_set)
    mt_splits = canonical_splits(treefiles["mitochondria"], shared_set)
    conflict_path = supplemental_dir / "strongly_supported_organelle_conflicts.tsv"
    conflict_state_path = state_dir / "conflicts.json"
    conflict_fingerprint = build_stage_fingerprint_from_hashes(
        "tree:organelle_conflicts",
        {
            **runtime_provenance(
                root,
                {
                    "biopython": ("python", "-c", "import Bio; print(Bio.__version__)"),
                },
            ),
            config_path.relative_to(root).as_posix(): sha256_file(config_path),
            treefiles["chloroplast"].relative_to(root).as_posix(): sha256_file(treefiles["chloroplast"]),
            treefiles["mitochondria"].relative_to(root).as_posix(): sha256_file(treefiles["mitochondria"]),
            (metadata_dir / "shared_samples.tsv").relative_to(root).as_posix(): sha256_file(metadata_dir / "shared_samples.tsv"),
        },
        {
            "chloroplast": json.loads((state_dir / "chloroplast.json").read_text())["fingerprint"]["digest"],
            "mitochondria": json.loads((state_dir / "mitochondria.json").read_text())["fingerprint"]["digest"],
        },
        [
            "report incompatible unrooted splits only when "
            f"SH-aLRT>={minimum_strong_sh_alrt:g} and UFBoot>={minimum_strong_ufboot:g} in both organelles"
        ],
    )
    if args.resume and conflict_state_path.exists():
        saved_conflicts = json.loads(conflict_state_path.read_text())
        validate_resume(saved_conflicts["fingerprint"]["digest"], conflict_fingerprint)
        validate_saved_outputs(root, saved_conflicts)
    else:
        if conflict_state_path.exists() or conflict_path.exists():
            raise RuntimeError("Existing unvalidated organelle-conflict output")
        conflict_count = 0
        with conflict_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "chloroplast_split",
                    "chloroplast_sh_alrt",
                    "chloroplast_ultrafast_bootstrap",
                    "mitochondria_split",
                    "mitochondria_sh_alrt",
                    "mitochondria_ultrafast_bootstrap",
                ]
            )
            for cp_split, cp_support in sorted(cp_splits.items(), key=lambda item: sorted(item[0])):
                if not is_strong_iqtree_support(cp_support, minimum_strong_sh_alrt, minimum_strong_ufboot):
                    continue
                for mt_split, mt_support in sorted(mt_splits.items(), key=lambda item: sorted(item[0])):
                    if is_strong_iqtree_support(
                        mt_support,
                        minimum_strong_sh_alrt,
                        minimum_strong_ufboot,
                    ) and incompatible(cp_split, mt_split, shared_set):
                        writer.writerow(
                            [
                                ",".join(sorted(cp_split)),
                                cp_support.sh_alrt,
                                cp_support.ultrafast_bootstrap,
                                ",".join(sorted(mt_split)),
                                mt_support.sh_alrt,
                                mt_support.ultrafast_bootstrap,
                            ]
                        )
                        conflict_count += 1
        conflict_state_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "strong_conflict_count": conflict_count,
                    "fingerprint": asdict(conflict_fingerprint),
                    "outputs": {conflict_path.relative_to(root).as_posix(): sha256_file(conflict_path)},
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
    print("completed supplementary partitioned concatenation and conflict report")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
