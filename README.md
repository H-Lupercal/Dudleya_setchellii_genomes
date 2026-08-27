# Dudleya organelle population genomics

This repository is organized by scientific status:

- [`source_data/`](source_data/) — immutable raw reads and reference candidates.
- [`canonical_publication/`](canonical_publication/) — repaired code and accepted publication outputs.
- [`canonical_publication/provenance/archive/`](canonical_publication/provenance/archive/) — audit metadata for pre-remediation artifacts; tracked snapshot contents are preserved on `archive/noncanonical-2026-08-17`.

Accepted scientific run: `publication-20260817` (acceptance PASS). Its historical
acceptance record, artifact manifest, and `CURRENT_RUN` pointer remain unchanged.

Current publication package: `publication-20260817-package-20260826`
(`PUBLICATION_PACKAGE_PASS`). It re-attests the current repository packaging
against the unchanged accepted science; it does not claim that the analysis was
rerun. See `canonical_publication/CURRENT_PACKAGE` and verify it from any clean
`main` checkout with:

```bash
PYTHONPATH=canonical_publication/pipeline/src \
  python canonical_publication/pipeline/scripts/attest_publication_package.py \
  verify --repository-root .
```

The sole supported scientific entrypoint remains
`canonical_publication/run_pipeline.sh`.
