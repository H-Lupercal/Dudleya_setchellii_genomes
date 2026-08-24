"""Dependency-ordered supplementary workflow."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tomllib
from pathlib import Path

from .configuration import validate_config

STAGES = (
    "canonical_guard",
    "metadata",
    "identity",
    "sensitivity",
    "claims",
    "inheritance",
    "phase1_gate",
    "likelihood_mapping",
    "comparative_analyses",
    "figures",
    "reports",
    "acceptance",
    "canonical_guard_final",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--from-stage", choices=STAGES, default=STAGES[0])
    parser.add_argument("--until-stage", choices=STAGES, default=STAGES[-1])
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parents[4]
    config_path = (root / args.config).resolve() if not args.config.is_absolute() else args.config.resolve()
    if not config_path.is_relative_to(root / "supplementary_analysis"):
        raise SystemExit("Supplementary configuration must be beneath supplementary_analysis/")
    config = tomllib.loads(config_path.read_text())
    validate_config(config)
    start, end = STAGES.index(args.from_stage), STAGES.index(args.until_stage)
    if start > end:
        raise SystemExit("--from-stage must not follow --until-stage")
    if start and not args.resume:
        raise SystemExit("--from-stage requires --resume")
    runner = root / "supplementary_analysis/pipeline/scripts/run_stage.py"
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    environment["PYTHONPYCACHEPREFIX"] = str(root / "supplementary_analysis/work/pycache")
    environment["MPLCONFIGDIR"] = str(root / "supplementary_analysis/work/matplotlib")
    supplement_source = str(root / "supplementary_analysis/pipeline/src")
    canonical_source = str(root / "canonical_publication/pipeline/src")
    environment["PYTHONPATH"] = os.pathsep.join([supplement_source, canonical_source, environment.get("PYTHONPATH", "")]).rstrip(os.pathsep)
    environment["PATH"] = str(root / ".tools/bioconda-env/bin") + os.pathsep + environment["PATH"]
    for stage in STAGES[: end + 1]:
        role = "UPSTREAM VALIDATION" if STAGES.index(stage) < start else "STAGE"
        command = [
            sys.executable,
            str(runner),
            "--stage",
            stage,
            "--config",
            str(config_path),
            "--run-id",
            args.run_id,
        ]
        if args.resume:
            command.append("--resume")
        print(f"{role} {stage}: {' '.join(command)}", flush=True)
        if not args.dry_run and STAGES.index(stage) >= start:
            subprocess.run(command, cwd=root, env=environment, check=True)
    return 0
