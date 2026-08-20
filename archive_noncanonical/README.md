# Noncanonical archive

Everything in the dated snapshot predates the publication remediation begun
on 2026-08-17. It is preserved for auditability and historical comparison
only. It must not be used as input to the canonical analysis.

The byte-preserving snapshot is in
`2026-08-17_pre_remediation/snapshot/`. `manifest.tsv` records each file or
symbolic link, its original path, size, SHA-256 checksum, tracking status, and
the reason for quarantine.

Canonical code rejects resolved paths beneath this directory, including
symbolic links that point here.

`manifest_tools/` contains archive-maintenance code and is not part of the
canonical scientific pipeline. The dated manifest inventories only the
byte-preserved `snapshot/` tree.
