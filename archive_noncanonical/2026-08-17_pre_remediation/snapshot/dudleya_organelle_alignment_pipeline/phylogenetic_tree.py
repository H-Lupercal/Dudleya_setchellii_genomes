"""Build cpDNA and mtDNA phylogenetic trees from callable consensus alignments.

This stage consumes the full callable-site FASTA
alignments and runs IQ-TREE maximum-likelihood tree inference separately for
cpDNA and mtDNA.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from dudleya_organelle_alignment_pipeline.pilot_alignment import shlex_join
from dudleya_organelle_alignment_pipeline.variant_calling import labeled_output_name


DEFAULT_CONSENSUS_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/11_callable_consensus"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/12_phylogenetic_tree"
)
DEFAULT_RUN_LABEL = "primary"
DEFAULT_MODEL = "GTR+F+G4"
DEFAULT_THREADS = 4
DEFAULT_FAST_SEARCH = True
DEFAULT_BOOTSTRAP_REPLICATES = 0


class PhylogeneticTreeError(RuntimeError):
    """Raised when this stage cannot safely build phylogenetic trees."""


@dataclass(frozen=True)
class TreeInput:
    organelle: str
    track_id: str
    sample_count: int
    alignment_sites: int
    missing_bases: int
    alignment_fasta_path: Path

    def to_result(
        self,
        model: str,
        method: str,
        tree_prefix: Path,
        treefile_path: Path,
        log_path: Path,
        iqtree_report_path: Path,
    ) -> "TreeResult":
        return TreeResult(
            organelle=self.organelle,
            track_id=self.track_id,
            sample_count=self.sample_count,
            alignment_sites=self.alignment_sites,
            missing_bases=self.missing_bases,
            model=model,
            method=method,
            alignment_fasta_path=self.alignment_fasta_path,
            tree_prefix=tree_prefix,
            treefile_path=treefile_path,
            log_path=log_path,
            iqtree_report_path=iqtree_report_path,
        )


@dataclass(frozen=True)
class TreeResult:
    organelle: str
    track_id: str
    sample_count: int
    alignment_sites: int
    missing_bases: int
    model: str
    method: str
    alignment_fasta_path: Path
    tree_prefix: Path
    treefile_path: Path
    log_path: Path
    iqtree_report_path: Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def require_iqtree() -> str:
    executable = shutil.which("iqtree") or shutil.which("iqtree2")
    if executable is None:
        raise PhylogeneticTreeError(
            "Missing required tool: iqtree or iqtree2. Activate the pipeline environment first."
        )
    return executable


def read_tree_inputs(consensus_dir: Path, run_label: str) -> list[TreeInput]:
    summary_path = consensus_dir / labeled_output_name(
        "callable_consensus_summary.tsv",
        run_label,
    )
    inputs: list[TreeInput] = []
    for row in read_tsv(summary_path):
        fasta_path = Path(row["alignment_fasta_path"])
        if not fasta_path.exists():
            raise PhylogeneticTreeError(f"Missing callable consensus FASTA: {fasta_path}")
        inputs.append(
            TreeInput(
                organelle=row["organelle"],
                track_id=row["track_id"],
                sample_count=int(row["sample_count"]),
                alignment_sites=int(row["consensus_length"]),
                missing_bases=int(row["missing_bases"]),
                alignment_fasta_path=fasta_path,
            )
        )
    if not inputs:
        raise PhylogeneticTreeError(f"No consensus rows found in {summary_path}")
    return inputs


def tree_output_prefix(tree_input: TreeInput, output_dir: Path, run_label: str) -> Path:
    label = f".{run_label}" if run_label else ""
    return output_dir / f"{tree_input.organelle}{label}.iqtree_ml"


def build_iqtree_command(
    iqtree_executable: str,
    alignment_path: Path,
    prefix: Path,
    model: str,
    threads: int,
    fast: bool,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
) -> list[str]:
    command = [
        iqtree_executable,
        "-s",
        alignment_path.as_posix(),
        "--seqtype",
        "DNA",
        "-m",
        model,
        "--prefix",
        prefix.as_posix(),
        "-T",
        str(threads),
        "--safe",
        "--redo",
        "--quiet",
    ]
    if fast:
        command.append("--fast")
    if bootstrap_replicates:
        command.extend(["-B", str(bootstrap_replicates), "--bnni"])
    return command


def tree_outputs_ready(prefix: Path) -> bool:
    treefile = Path(f"{prefix}.treefile")
    return treefile.exists() and treefile.stat().st_size > 0


def run_one_tree(
    tree_input: TreeInput,
    output_dir: Path,
    run_label: str,
    iqtree_executable: str,
    model: str,
    threads: int,
    fast: bool,
    bootstrap_replicates: int,
    force: bool,
) -> tuple[TreeResult, list[dict[str, str]]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = tree_output_prefix(tree_input, output_dir, run_label)
    command = build_iqtree_command(
        iqtree_executable=iqtree_executable,
        alignment_path=tree_input.alignment_fasta_path,
        prefix=prefix,
        model=model,
        threads=threads,
        fast=fast,
        bootstrap_replicates=bootstrap_replicates,
    )
    if bootstrap_replicates:
        method = f"iqtree_ml_ufboot{bootstrap_replicates}"
    else:
        method = "iqtree_ml_fast" if fast else "iqtree_ml"
    command_rows = [
        {
            "organelle": tree_input.organelle,
            "method": method,
            "command": shlex_join(command),
        }
    ]
    if tree_outputs_ready(prefix) and not force:
        command_rows.append(
            {
                "organelle": tree_input.organelle,
                "method": method,
                "command": "outputs already present; pass --force to regenerate",
            }
        )
    else:
        completed = subprocess.run(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if completed.returncode:
            raise PhylogeneticTreeError(
                f"IQ-TREE failed for {tree_input.organelle}; see {prefix}.log"
            )
    result = tree_input.to_result(
        model=model,
        method=method,
        tree_prefix=prefix,
        treefile_path=Path(f"{prefix}.treefile"),
        log_path=Path(f"{prefix}.log"),
        iqtree_report_path=Path(f"{prefix}.iqtree"),
    )
    if not result.treefile_path.exists() or result.treefile_path.stat().st_size == 0:
        raise PhylogeneticTreeError(
            f"Missing IQ-TREE treefile for {tree_input.organelle}: {result.treefile_path}"
        )
    return result, command_rows


def run_phylogenetic_trees(
    consensus_dir: Path = DEFAULT_CONSENSUS_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
    model: str = DEFAULT_MODEL,
    threads: int = DEFAULT_THREADS,
    fast: bool = DEFAULT_FAST_SEARCH,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    force: bool = False,
) -> list[TreeResult]:
    iqtree_executable = require_iqtree()
    inputs = read_tree_inputs(consensus_dir, run_label)
    results: list[TreeResult] = []
    command_rows: list[dict[str, str]] = []
    for tree_input in inputs:
        result, rows = run_one_tree(
            tree_input=tree_input,
            output_dir=output_dir,
            run_label=run_label,
            iqtree_executable=iqtree_executable,
            model=model,
            threads=threads,
            fast=fast,
            bootstrap_replicates=bootstrap_replicates,
            force=force,
        )
        results.append(result)
        command_rows.extend(rows)
    write_tree_outputs(output_dir, results, command_rows, run_label)
    return results


def write_tree_outputs(
    output_dir: Path,
    results: list[TreeResult],
    command_rows: list[dict[str, str]],
    run_label: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / labeled_output_name("phylogenetic_tree_summary.tsv", run_label),
        [
            {
                "organelle": result.organelle,
                "track_id": result.track_id,
                "sample_count": str(result.sample_count),
                "alignment_sites": str(result.alignment_sites),
                "missing_bases": str(result.missing_bases),
                "method": result.method,
                "model": result.model,
                "alignment_fasta_path": result.alignment_fasta_path.as_posix(),
                "treefile_path": result.treefile_path.as_posix(),
                "log_path": result.log_path.as_posix(),
                "iqtree_report_path": result.iqtree_report_path.as_posix(),
            }
            for result in results
        ],
        [
            "organelle",
            "track_id",
            "sample_count",
            "alignment_sites",
            "missing_bases",
            "method",
            "model",
            "alignment_fasta_path",
            "treefile_path",
            "log_path",
            "iqtree_report_path",
        ],
    )
    write_tsv(
        output_dir / labeled_output_name("phylogenetic_tree_commands.tsv", run_label),
        command_rows,
        ["organelle", "method", "command"],
    )
    write_report(
        output_dir / labeled_output_name("phylogenetic_tree_report.md", run_label),
        results=results,
        run_label=run_label,
    )


def write_report(path: Path, results: list[TreeResult], run_label: str) -> None:
    label = run_label or "full"
    lines = [
        "# Phylogenetic Trees",
        "",
        "This step builds cpDNA and mtDNA phylogenetic trees from the",
        "full callable-site consensus FASTA alignments using IQ-TREE maximum-likelihood",
        "inference. Bootstrap support is included when",
        "requested for the run.",
        "",
        "## Run",
        "",
        f"- Run label: `{label}`",
        f"- Method: {results[0].method if results else 'IQ-TREE maximum-likelihood'}",
        "",
        "## Results",
        "",
    ]
    for result in results:
        lines.extend(
            [
                f"### {result.organelle}",
                "",
                f"- Track: `{result.track_id}`",
                f"- Samples: {result.sample_count}",
                f"- Alignment sites: {result.alignment_sites}",
                f"- Missing bases: {result.missing_bases}",
                f"- Model: `{result.model}`",
                f"- Tree: `{result.treefile_path}`",
                f"- IQ-TREE report: `{result.iqtree_report_path}`",
                f"- Log: `{result.log_path}`",
                "",
            ]
        )
    path.write_text("\n".join(lines))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build cpDNA/mtDNA phylogenetic trees with IQ-TREE."
    )
    parser.add_argument("--consensus-dir", type=Path, default=DEFAULT_CONSENSUS_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--full-search", action="store_true")
    parser.add_argument("--bootstrap-replicates", type=int, default=DEFAULT_BOOTSTRAP_REPLICATES)
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    results = run_phylogenetic_trees(
        consensus_dir=args.consensus_dir,
        output_dir=args.output_dir,
        run_label=args.run_label,
        model=args.model,
        threads=args.threads,
        fast=not args.full_search and not args.bootstrap_replicates,
        bootstrap_replicates=args.bootstrap_replicates,
        force=args.force,
    )
    for result in results:
        print(
            f"{result.organelle}: {result.method} tree for "
            f"{result.sample_count} samples at {result.treefile_path}"
        )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
