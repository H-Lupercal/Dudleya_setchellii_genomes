#!/usr/bin/env python3
"""Run Step 7: call raw haploid cpDNA/mtDNA variants."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.variant_calling import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
