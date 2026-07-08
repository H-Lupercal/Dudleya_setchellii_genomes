#!/usr/bin/env python3
"""Run: cpDNA/mtDNA ADMIXTURE-style clustering."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.admixture_analysis import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
