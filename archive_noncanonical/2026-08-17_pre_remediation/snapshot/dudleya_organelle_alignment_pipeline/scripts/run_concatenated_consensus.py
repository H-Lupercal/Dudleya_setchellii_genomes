#!/usr/bin/env python3
"""Run: concatenate matched cpDNA and mtDNA consensus alignments."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.concatenated_consensus import (  # noqa: E402
    main,
)


if __name__ == "__main__":
    raise SystemExit(main())
