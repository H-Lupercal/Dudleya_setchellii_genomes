#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from dudleya_supplement.stages import run_stage
from dudleya_supplement.workflow import STAGES


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=STAGES, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[3]
    run_stage(args.stage, root, args.config.resolve(), args.run_id, args.resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
