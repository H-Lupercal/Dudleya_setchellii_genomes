# Dudleya organelle population genomics

This repository is organized by scientific status:

- [`source_data/`](source_data/) — immutable raw reads and reference candidates.
- [`canonical_publication/`](canonical_publication/) — repaired code and accepted publication outputs.
- [`canonical_publication/provenance/archive/`](canonical_publication/provenance/archive/) — audit metadata for pre-remediation artifacts; tracked snapshot contents are preserved on `archive/noncanonical-2026-08-17`.

Current canonical run: `publication-20260817` (acceptance PASS). The sole supported entrypoint is `canonical_publication/run_pipeline.sh`; see `canonical_publication/CURRENT_RUN` and `canonical_publication/provenance/runs/` for checksummed status.
