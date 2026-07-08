#!/usr/bin/env python3
"""Run step 4: create cpDNA/mtDNA analysis masks and track manifests."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.analysis_masks import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
