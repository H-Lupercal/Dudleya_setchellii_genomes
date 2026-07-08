#!/usr/bin/env python3
"""Run Step 9: build SNP-only cpDNA/mtDNA FASTA alignments."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.snp_alignment import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
