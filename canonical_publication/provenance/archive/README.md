# Noncanonical archive provenance

The tracked pre-remediation snapshot is preserved on the remote branch
`archive/noncanonical-2026-08-17`. It must not be used as input to the
canonical scientific analysis.

Preservation identifiers:

- Source commit: `abb16527d20a9dd949261d5ca2bc602987a82cee`
- `archive_noncanonical` subtree: `6d2f2ed95021de132c561017008b78cb47a3a294`
- Tracked files beneath that subtree: 1,717
- Manifest SHA-256: `7d7d0eb52daaf27c0d12f0608d37b064e15c8161fc5d79b86cd19a828a7ef047`

The byte-identical manifest is retained at
`2026-08-17_pre_remediation/manifest.tsv` so canonical reporting can account
for all 5,674 historical entries without checking out the snapshot. The
manifest also inventories workstation-only files that were ignored by Git;
those approximately 69 GB of local outputs are not backed up by the remote
branch. A full filesystem checksum validation therefore requires the original
local-only files in addition to the archive branch.

Existing accepted-run provenance files are unchanged historical records. Paths
within them are relative to the source commit above, where the original archive
layout remains available.
