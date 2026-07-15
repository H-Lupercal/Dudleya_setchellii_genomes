#!/usr/bin/env python3
"""Run Stage 21: build cpDNA and mtDNA haplotype networks."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.haplotype_network import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
