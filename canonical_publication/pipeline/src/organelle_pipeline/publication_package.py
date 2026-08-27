"""Deterministic attestation for the current canonical publication package."""

from __future__ import annotations

import csv
import hashlib
import io
import json
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from .paths import validate_run_id
from .provenance import pipeline_code_digest, sha256_file

PACKAGE_MANIFEST_FIELDS = ("canonical_path", "size_bytes", "sha256", "status")
PACKAGE_STATUSES = frozenset({"canonical", "archive_audit_manifest", "historical_provenance"})
ARCHIVE_MANIFEST_OLD_PATH = "archive_noncanonical/2026-08-17_pre_remediation/manifest.tsv"
ARCHIVE_MANIFEST_PATH = "canonical_publication/provenance/archive/2026-08-17_pre_remediation/manifest.tsv"
ARCHIVE_SNAPSHOT = {
    "branch": "archive/noncanonical-2026-08-17",
    "source_commit": "abb16527d20a9dd949261d5ca2bc602987a82cee",
    "subtree": "6d2f2ed95021de132c561017008b78cb47a3a294",
    "tracked_file_count": 1717,
    "manifest_sha256": "7d7d0eb52daaf27c0d12f0608d37b064e15c8161fc5d79b86cd19a828a7ef047",
    "scope": "tracked files only; ignored local outputs excluded",
}
KNOWN_BASE_EVIDENCE = {
    "publication-20260817": {
        "acceptance_sha256": "585f72eef82130f5709aed9a850367c76dd8035945a18ddc5d86e3ab1c3c4fda",
        "artifact_manifest_sha256": "0817889739b03b9b33c4fe3aa41baf95d6d6b20b2a2a2daf314a7dbf952890f8",
    }
}
MUTABLE_HISTORICAL_PATHS = frozenset({"README.md", "canonical_publication/README.md"})
MUTABLE_HISTORICAL_PREFIXES = ("canonical_publication/pipeline/",)
IGNORED_DIRECTORY_NAMES = frozenset({"__pycache__", ".pytest_cache", ".ruff_cache", ".matplotlib"})
ACCEPTANCE_FIELDS = frozenset(
    {
        "acceptance_scope",
        "archive_snapshot",
        "artifact_count",
        "artifact_manifest",
        "artifact_manifest_sha256",
        "base_run",
        "errors",
        "package_id",
        "pipeline_code_digest",
        "status",
        "status_label",
    }
)
BASE_EVIDENCE_FIELDS = frozenset(
    {
        "run_id",
        "acceptance_path",
        "acceptance_sha256",
        "artifact_manifest_path",
        "artifact_manifest_sha256",
        "historical_artifact_count",
    }
)


class PublicationPackageError(RuntimeError):
    """Raised when package creation or verification cannot establish integrity."""


@dataclass(frozen=True)
class PackagePaths:
    manifest: Path
    acceptance: Path
    current: Path
    checksum_index: Path


def _package_paths(root: Path, package_id: str) -> PackagePaths:
    manifest_dir = root / "canonical_publication/provenance/manifests"
    return PackagePaths(
        manifest=manifest_dir / f"{package_id}.final_artifacts.tsv",
        acceptance=root / f"canonical_publication/provenance/packages/{package_id}/ACCEPTANCE.json",
        current=root / "canonical_publication/CURRENT_PACKAGE",
        checksum_index=manifest_dir / f"{package_id}.acceptance.sha256",
    )


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.is_symlink():
        raise PublicationPackageError(f"Required TSV is missing or not a regular file: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != PACKAGE_MANIFEST_FIELDS:
            raise PublicationPackageError(f"Package manifest schema mismatch: {path}")
        return list(reader)


def _validated_metadata(row: dict[str, str], row_number: int, seen_paths: set[str]) -> tuple[str, int, str, str]:
    if tuple(row) != PACKAGE_MANIFEST_FIELDS:
        raise PublicationPackageError(f"Package manifest schema mismatch in row {row_number}")
    relative = row["canonical_path"]
    pure_path = PurePosixPath(relative)
    if not relative or relative == "." or pure_path.is_absolute() or ".." in pure_path.parts or "\\" in relative:
        raise PublicationPackageError(f"Unsafe package path: {relative}")
    if relative in seen_paths:
        raise PublicationPackageError(f"Duplicate package path: {relative}")
    seen_paths.add(relative)
    size_text = row["size_bytes"]
    try:
        size = int(size_text)
    except ValueError as error:
        raise PublicationPackageError(f"Invalid package size: {size_text}") from error
    if size < 0 or str(size) != size_text:
        raise PublicationPackageError(f"Invalid package size: {size_text}")
    digest = row["sha256"]
    if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
        raise PublicationPackageError(f"Invalid package SHA-256: {digest}")
    status = row["status"]
    if status not in PACKAGE_STATUSES:
        raise PublicationPackageError(f"Invalid package status: {status}")
    return relative, size, digest, status


def _verify_regular_file(root: Path, relative: str, size: int, digest: str) -> None:
    pure_path = PurePosixPath(relative)
    manifest_path = root / Path(*pure_path.parts)
    resolved = manifest_path.resolve(strict=False)
    if resolved != root and not resolved.is_relative_to(root):
        raise PublicationPackageError(f"Package path escapes repository: {relative}")
    if manifest_path.is_symlink() or not resolved.is_file():
        raise PublicationPackageError(f"Package path is not a regular file: {relative}")
    if resolved.stat().st_size != size:
        raise PublicationPackageError(f"Package size mismatch: {relative}")
    if sha256_file(resolved) != digest:
        raise PublicationPackageError(f"Package checksum mismatch: {relative}")


def _base_run_evidence(root: Path, base_run_id: str) -> dict[str, str | int]:
    validate_run_id(base_run_id)
    acceptance_path = root / f"canonical_publication/provenance/runs/{base_run_id}/ACCEPTANCE.json"
    manifest_path = root / f"canonical_publication/provenance/manifests/{base_run_id}.final_artifacts.tsv"
    try:
        acceptance = json.loads(acceptance_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise PublicationPackageError(f"Base-run acceptance is unreadable: {acceptance_path}") from error
    expected_manifest = manifest_path.relative_to(root).as_posix()
    if acceptance.get("status") != "PASS" or acceptance.get("run_id") != base_run_id:
        raise PublicationPackageError(f"Base run is not accepted: {base_run_id}")
    if acceptance.get("final_manifest") != expected_manifest:
        raise PublicationPackageError(f"Base-run manifest pointer mismatch: {base_run_id}")

    acceptance_digest = sha256_file(acceptance_path)
    manifest_digest = sha256_file(manifest_path)
    known = KNOWN_BASE_EVIDENCE.get(base_run_id)
    if known and acceptance_digest != known["acceptance_sha256"]:
        raise PublicationPackageError(f"Historical acceptance checksum mismatch: {base_run_id}")
    if known and manifest_digest != known["artifact_manifest_sha256"]:
        raise PublicationPackageError(f"Historical artifact manifest checksum mismatch: {base_run_id}")

    rows = _read_tsv(manifest_path)
    seen_paths: set[str] = set()
    for row_number, row in enumerate(rows, start=1):
        relative, size, digest, _ = _validated_metadata(row, row_number, seen_paths)
        if relative == ARCHIVE_MANIFEST_OLD_PATH:
            _verify_regular_file(root, ARCHIVE_MANIFEST_PATH, size, digest)
            if known and digest != ARCHIVE_SNAPSHOT["manifest_sha256"]:
                raise PublicationPackageError("Historical archive manifest checksum is not the recorded snapshot")
            continue
        if relative in MUTABLE_HISTORICAL_PATHS or relative.startswith(MUTABLE_HISTORICAL_PREFIXES):
            mutable_path = root / relative
            if mutable_path.is_symlink() or not mutable_path.is_file():
                raise PublicationPackageError(f"Missing historical packaging file: {relative}")
            continue
        _verify_regular_file(root, relative, size, digest)
    return {
        "run_id": base_run_id,
        "acceptance_path": acceptance_path.relative_to(root).as_posix(),
        "acceptance_sha256": acceptance_digest,
        "artifact_manifest_path": expected_manifest,
        "artifact_manifest_sha256": manifest_digest,
        "historical_artifact_count": len(rows),
    }


def _is_excluded(path: Path, root: Path, outputs: PackagePaths) -> bool:
    if path in {outputs.manifest, outputs.acceptance, outputs.current, outputs.checksum_index}:
        return True
    relative = path.relative_to(root)
    if "work" in relative.parts and relative.parts[:2] == ("canonical_publication", "work"):
        return True
    return bool(IGNORED_DIRECTORY_NAMES & set(relative.parts)) or path.suffix == ".pyc"


def _artifact_status(relative: str) -> str:
    if relative.startswith("canonical_publication/provenance/archive/"):
        return "archive_audit_manifest"
    if relative.startswith("canonical_publication/provenance/"):
        return "historical_provenance"
    return "canonical"


def _current_artifact_rows(root: Path, outputs: PackagePaths, current_content: bytes) -> list[dict[str, str]]:
    candidates: set[Path] = set()
    canonical = root / "canonical_publication"
    if not canonical.is_dir():
        raise PublicationPackageError("canonical_publication is missing")
    for path in canonical.rglob("*"):
        if path.is_symlink() and not _is_excluded(path, root, outputs):
            raise PublicationPackageError(f"Publication package cannot contain a symlink: {path.relative_to(root)}")
        if path.is_file() and not _is_excluded(path, root, outputs):
            candidates.add(path)
    for relative in ("README.md", ".gitignore", ".gitattributes"):
        path = root / relative
        if path.is_file() and not path.is_symlink():
            candidates.add(path)
    workflows = root / ".github/workflows"
    if workflows.is_dir():
        candidates.update(path for path in workflows.rglob("*") if path.is_file() and not path.is_symlink())

    rows = []
    for path in candidates:
        relative = path.relative_to(root).as_posix()
        rows.append(
            {
                "canonical_path": relative,
                "size_bytes": str(path.stat().st_size),
                "sha256": sha256_file(path),
                "status": _artifact_status(relative),
            }
        )
    current_relative = outputs.current.relative_to(root).as_posix()
    rows.append(
        {
            "canonical_path": current_relative,
            "size_bytes": str(len(current_content)),
            "sha256": hashlib.sha256(current_content).hexdigest(),
            "status": "canonical",
        }
    )
    return sorted(rows, key=lambda row: row["canonical_path"])


def _tsv_bytes(rows: list[dict[str, str]]) -> bytes:
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=PACKAGE_MANIFEST_FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return handle.getvalue().encode()


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def _write_new_or_identical(path: Path, content: bytes) -> None:
    if path.exists():
        if path.is_symlink() or not path.is_file() or path.read_bytes() != content:
            raise PublicationPackageError(f"Refusing to replace differing package artifact: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def create_publication_package(
    repository_root: Path | str,
    base_run_id: str,
    package_id: str,
) -> PackagePaths:
    """Create a deterministic package attestation."""

    root = Path(repository_root).resolve()
    validate_run_id(package_id)
    outputs = _package_paths(root, package_id)
    base_evidence = _base_run_evidence(root, base_run_id)
    current_content = f"{package_id}\tPASS\n".encode()
    rows = _current_artifact_rows(root, outputs, current_content)
    manifest_content = _tsv_bytes(rows)
    manifest_digest = hashlib.sha256(manifest_content).hexdigest()
    archive_snapshot = dict(ARCHIVE_SNAPSHOT)
    archive_copy = root / ARCHIVE_MANIFEST_PATH
    if not archive_copy.is_file() or archive_copy.is_symlink():
        raise PublicationPackageError("Canonical archive manifest copy is missing")
    archive_snapshot["manifest_sha256"] = sha256_file(archive_copy)
    if base_run_id in KNOWN_BASE_EVIDENCE and archive_snapshot["manifest_sha256"] != ARCHIVE_SNAPSHOT["manifest_sha256"]:
        raise PublicationPackageError("Canonical archive manifest checksum mismatch")
    acceptance = {
        "acceptance_scope": "publication_package_and_provenance",
        "archive_snapshot": archive_snapshot,
        "artifact_count": len(rows),
        "artifact_manifest": outputs.manifest.relative_to(root).as_posix(),
        "artifact_manifest_sha256": manifest_digest,
        "base_run": base_evidence,
        "errors": [],
        "package_id": package_id,
        "pipeline_code_digest": pipeline_code_digest(root),
        "status": "PASS",
        "status_label": "PUBLICATION_PACKAGE_PASS",
    }
    acceptance_content = _json_bytes(acceptance)
    acceptance_digest = hashlib.sha256(acceptance_content).hexdigest()
    checksum_content = (
        f"{manifest_digest}  {outputs.manifest.relative_to(root).as_posix()}\n"
        f"{acceptance_digest}  {outputs.acceptance.relative_to(root).as_posix()}\n"
        f"{hashlib.sha256(current_content).hexdigest()}  {outputs.current.relative_to(root).as_posix()}\n"
    ).encode()
    desired = {
        outputs.current: current_content,
        outputs.manifest: manifest_content,
        outputs.acceptance: acceptance_content,
        outputs.checksum_index: checksum_content,
    }
    differing = [path for path, content in desired.items() if path.exists() and (path.is_symlink() or path.read_bytes() != content)]
    if differing:
        raise PublicationPackageError(f"Refusing to replace differing package artifacts: {differing}")
    for path, content in desired.items():
        _write_new_or_identical(path, content)
    verify_publication_package(root, package_id)
    return outputs


def verify_publication_package(repository_root: Path | str, package_id: str | None = None) -> dict[str, object]:
    """Verify a previously created publication package."""

    root = Path(repository_root).resolve()
    current_path = root / "canonical_publication/CURRENT_PACKAGE"
    try:
        current_content = current_path.read_text()
    except OSError as error:
        raise PublicationPackageError("CURRENT_PACKAGE is missing") from error
    current_fields = current_content.removesuffix("\n").split("\t")
    if len(current_fields) != 2 or current_fields[1] != "PASS" or current_content != f"{current_fields[0]}\tPASS\n":
        raise PublicationPackageError("CURRENT_PACKAGE is malformed")
    current_package_id = current_fields[0]
    validate_run_id(current_package_id)
    if package_id is not None and package_id != current_package_id:
        raise PublicationPackageError("CURRENT_PACKAGE does not match the requested package")
    outputs = _package_paths(root, current_package_id)
    try:
        acceptance = json.loads(outputs.acceptance.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise PublicationPackageError("Package acceptance is unreadable") from error
    if not isinstance(acceptance, dict) or acceptance.keys() != ACCEPTANCE_FIELDS:
        raise PublicationPackageError("Package acceptance schema mismatch")
    if (
        acceptance.get("package_id") != current_package_id
        or acceptance.get("status") != "PASS"
        or acceptance.get("status_label") != "PUBLICATION_PACKAGE_PASS"
        or acceptance.get("errors") != []
    ):
        raise PublicationPackageError("Package acceptance does not record PASS")
    if acceptance.get("pipeline_code_digest") != pipeline_code_digest(root):
        raise PublicationPackageError("Package pipeline code digest mismatch")
    manifest_digest = sha256_file(outputs.manifest)
    if acceptance.get("artifact_manifest") != outputs.manifest.relative_to(root).as_posix():
        raise PublicationPackageError("Package artifact manifest pointer mismatch")
    if acceptance.get("artifact_manifest_sha256") != manifest_digest:
        raise PublicationPackageError("Package artifact manifest checksum mismatch")
    rows = _read_tsv(outputs.manifest)
    if acceptance.get("artifact_count") != validate_package_manifest(rows, root):
        raise PublicationPackageError("Package artifact count mismatch")
    base_run = acceptance.get("base_run")
    if not isinstance(base_run, dict) or base_run.keys() != BASE_EVIDENCE_FIELDS or not isinstance(base_run.get("run_id"), str):
        raise PublicationPackageError("Package base-run evidence is malformed")
    if _base_run_evidence(root, base_run["run_id"]) != base_run:
        raise PublicationPackageError("Package base-run evidence mismatch")
    archive_snapshot = acceptance.get("archive_snapshot")
    if not isinstance(archive_snapshot, dict) or archive_snapshot.keys() != ARCHIVE_SNAPSHOT.keys():
        raise PublicationPackageError("Package archive snapshot evidence is malformed")
    archive_copy_digest = sha256_file(root / ARCHIVE_MANIFEST_PATH)
    if archive_snapshot.get("manifest_sha256") != archive_copy_digest:
        raise PublicationPackageError("Package archive snapshot checksum mismatch")
    for key, expected in ARCHIVE_SNAPSHOT.items():
        if key != "manifest_sha256" and archive_snapshot.get(key) != expected:
            raise PublicationPackageError(f"Package archive snapshot {key} mismatch")

    expected_index = [
        (manifest_digest, outputs.manifest),
        (sha256_file(outputs.acceptance), outputs.acceptance),
        (sha256_file(outputs.current), outputs.current),
    ]
    try:
        index_lines = outputs.checksum_index.read_text().splitlines()
    except OSError as error:
        raise PublicationPackageError("Package checksum index is missing") from error
    if len(index_lines) != len(expected_index):
        raise PublicationPackageError("Package checksum index is malformed")
    for line, (expected_digest, expected_path) in zip(index_lines, expected_index, strict=True):
        if line != f"{expected_digest}  {expected_path.relative_to(root).as_posix()}":
            raise PublicationPackageError("Package checksum index is malformed or stale")
    return acceptance


def validate_package_manifest(rows: list[dict[str, str]], repository_root: Path | str) -> int:
    """Validate package rows against the repository and return the file count."""

    root = Path(repository_root).resolve()
    seen_paths: set[str] = set()
    for row_number, row in enumerate(rows, start=1):
        relative, size, digest, _ = _validated_metadata(row, row_number, seen_paths)
        _verify_regular_file(root, relative, size, digest)
    return len(rows)
