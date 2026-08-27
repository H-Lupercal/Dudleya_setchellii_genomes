#!/usr/bin/env python3
"""Create or verify a canonical publication-package attestation."""

from __future__ import annotations

import argparse
from pathlib import Path

from organelle_pipeline.publication_package import create_publication_package, verify_publication_package


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("create", help="create a new immutable package attestation")
    create.add_argument("--repository-root", type=Path, required=True)
    create.add_argument("--base-run-id", required=True)
    create.add_argument("--package-id", required=True)
    verify = commands.add_parser("verify", help="verify the current package without modifying files")
    verify.add_argument("--repository-root", type=Path, required=True)
    verify.add_argument("--package-id")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.command == "create":
        created = create_publication_package(args.repository_root, args.base_run_id, args.package_id)
        print(f"publication package created: {args.package_id} ({created.manifest})")
    else:
        acceptance = verify_publication_package(args.repository_root, args.package_id)
        print(f"publication package verified: {acceptance['package_id']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
