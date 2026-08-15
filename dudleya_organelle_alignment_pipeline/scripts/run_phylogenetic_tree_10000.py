#!/usr/bin/env python3
"""Run cpDNA and mtDNA trees with 10,000 bootstrap replicates each."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.phylogenetic_tree import (  # noqa: E402
    main as run_phylogenetic_tree_main,
)


OUTPUT_DIR = (
    "dudleya_organelle_alignment_pipeline/results/"
    "19_bootstrap_phylogenetic_tree_10000"
)


def build_run_arguments() -> list[str]:
    """Return the fixed settings for the 20,000-total-bootstrap run."""
    return [
        "--run-label",
        "primary",
        "--output-dir",
        OUTPUT_DIR,
        "--threads",
        "14",
        "--bootstrap-replicates",
        "10000",
    ]


def main() -> int:
    return run_phylogenetic_tree_main(build_run_arguments())


if __name__ == "__main__":
    raise SystemExit(main())
