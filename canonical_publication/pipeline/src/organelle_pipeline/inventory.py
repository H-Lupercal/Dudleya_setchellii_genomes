"""Checksummed filesystem inventories for archived and source artifacts."""

from __future__ import annotations

import csv
import hashlib
import os
from dataclasses import asdict, dataclass
from pathlib import Path

from .provenance import sha256_file

ACCEPTABLE_SOURCE_VALIDATION_STATUSES = frozenset(
    {
        "PASS",
        "PASS_WITH_DECLARED_MISSING",
        "PASS_WITH_PROVIDER_METADATA_WARNING",
        "PASS_WITH_DECLARED_MISSING_AND_PROVIDER_METADATA_WARNING",
    }
)


def classify_provider_md5(
    expected: str,
    observed: str,
    source_is_provider_manifest: bool,
    source_inventory_matches: bool,
) -> str:
    """Classify provider MD5 evidence without treating a manifest as self-authenticating."""

    if not source_inventory_matches:
        return "FAIL_CHECKSUM"
    if expected.lower() == observed.lower():
        return "PASS"
    if source_is_provider_manifest:
        return "UNVERIFIABLE_SELF_REFERENCE"
    return "FAIL_CHECKSUM"


def source_validation_status(
    has_failures: bool,
    has_declared_missing: bool,
    has_self_reference_warning: bool,
) -> str:
    """Summarize source validation without hiding provider-metadata warnings."""

    if has_failures:
        return "FAIL"
    if has_declared_missing and has_self_reference_warning:
        return "PASS_WITH_DECLARED_MISSING_AND_PROVIDER_METADATA_WARNING"
    if has_declared_missing:
        return "PASS_WITH_DECLARED_MISSING"
    if has_self_reference_warning:
        return "PASS_WITH_PROVIDER_METADATA_WARNING"
    return "PASS"


@dataclass(frozen=True)
class ArtifactRecord:
    original_path: str
    archived_path: str
    artifact_type: str
    size_bytes: int
    sha256: str
    git_status: str
    reason: str


def _symlink_digest(path: Path) -> str:
    return hashlib.sha256(os.readlink(path).encode()).hexdigest()


def inventory_tree(
    snapshot_root: Path | str,
    repository_root: Path | str,
    tracked_original_paths: set[str],
    reason: str,
) -> list[ArtifactRecord]:
    """Inventory files and links without following links outside the tree."""

    snapshot = Path(snapshot_root).resolve()
    repository = Path(repository_root).resolve()
    records: list[ArtifactRecord] = []
    for path in sorted(snapshot.rglob("*")):
        if path.is_symlink():
            artifact_type = "symlink"
            size = path.lstat().st_size
            digest = _symlink_digest(path)
        elif path.is_file():
            artifact_type = "file"
            size = path.stat().st_size
            digest = sha256_file(path)
        else:
            continue
        original = path.relative_to(snapshot).as_posix()
        archived = path.relative_to(repository).as_posix()
        records.append(
            ArtifactRecord(
                original_path=original,
                archived_path=archived,
                artifact_type=artifact_type,
                size_bytes=size,
                sha256=digest,
                git_status="tracked" if original in tracked_original_paths else "local-only",
                reason=reason,
            )
        )
    return records


def write_inventory(path: Path | str, records: list[ArtifactRecord]) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [field.name for field in ArtifactRecord.__dataclass_fields__.values()]
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(asdict(record) for record in records)


def validate_inventory(rows: list[dict[str, str]], repository_root: Path | str) -> str:
    """Verify every manifested file/link and return an aggregate content digest."""

    root = Path(repository_root).resolve()
    aggregate = hashlib.sha256()
    seen_paths: set[str] = set()
    for row in rows:
        relative = row["archived_path"]
        if relative in seen_paths:
            raise ValueError(f"Duplicate inventory path: {relative}")
        seen_paths.add(relative)
        relative_path = Path(relative)
        if relative_path.is_absolute() or ".." in relative_path.parts:
            raise ValueError(f"Inventory path escapes repository: {relative}")
        manifest_path = root / relative_path
        path = manifest_path.resolve(strict=False)
        expected_type = row["artifact_type"]
        if expected_type == "symlink":
            if not manifest_path.is_symlink():
                raise ValueError(f"Inventory type mismatch: {relative}")
            observed_size = manifest_path.lstat().st_size
            observed_digest = _symlink_digest(manifest_path)
        elif expected_type == "file":
            if path != root and not path.is_relative_to(root):
                raise ValueError(f"Inventory file resolves outside repository: {relative}")
            if not path.is_file() or (root / relative).is_symlink():
                raise ValueError(f"Inventory type mismatch: {relative}")
            observed_size = path.stat().st_size
            observed_digest = sha256_file(path)
        else:
            raise ValueError(f"Unsupported inventory artifact type: {expected_type}")
        if observed_size != int(row["size_bytes"]):
            raise ValueError(f"Inventory size mismatch: {relative}")
        if observed_digest != row["sha256"]:
            raise ValueError(f"Inventory checksum mismatch: {relative}")
        aggregate.update(f"{relative}\0{expected_type}\0{observed_size}\0{observed_digest}\n".encode())
    return aggregate.hexdigest()
