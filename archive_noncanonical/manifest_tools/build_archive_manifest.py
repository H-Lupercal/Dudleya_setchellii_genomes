#!/usr/bin/env python3
"""Maintenance tool for rebuilding the checksummed legacy-snapshot manifest."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from organelle_pipeline.inventory import inventory_tree, write_inventory


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    completed = subprocess.run(
        ["git", "ls-files"],
        cwd=args.repository_root,
        check=True,
        capture_output=True,
        text=True,
    )
    tracked_paths = set(completed.stdout.splitlines())
    snapshot_relative = (
        args.snapshot.resolve().relative_to(args.repository_root.resolve()).as_posix()
    )
    tracked = {
        path.removeprefix(f"{snapshot_relative}/")
        if path.startswith(f"{snapshot_relative}/")
        else path
        for path in tracked_paths
    }
    records = inventory_tree(
        args.snapshot,
        repository_root=args.repository_root,
        tracked_original_paths=tracked,
        reason="pre-remediation noncanonical artifact; preserved for audit only",
    )
    write_inventory(args.output, records)
    print(f"wrote {len(records)} artifact records to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
