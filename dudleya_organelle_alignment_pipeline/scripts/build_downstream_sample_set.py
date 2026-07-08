#!/usr/bin/env python3
"""Run: build downstream include/exclude sample sets."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.downstream_sample_set import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
