#!/usr/bin/env python3
"""Run step 5: all-sample cpDNA/mtDNA alignment and track-aware QC."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.all_sample_alignment import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
