#!/usr/bin/env python3
"""Build SHA-256 inventory for immutable source inputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from organelle_pipeline.inventory import inventory_tree
from organelle_pipeline.paths import repository_relative


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.repository_root.resolve()
    source_root = root / repository_relative(args.source_root, root)
    output = root / repository_relative(args.output, root)
    records = inventory_tree(
        source_root,
        repository_root=root,
        tracked_original_paths=set(),
        reason="immutable source input",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["source_path", "artifact_type", "size_bytes", "sha256", "immutable"])
        for record in records:
            writer.writerow(
                [
                    record.archived_path,
                    record.artifact_type,
                    record.size_bytes,
                    record.sha256,
                    "yes",
                ]
            )
    print(f"wrote {len(records)} source records to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
