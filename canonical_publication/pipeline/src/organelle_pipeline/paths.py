"""Repository path guards for canonical analysis code."""

from __future__ import annotations

import re
from pathlib import Path


class CanonicalPathError(ValueError):
    """Raised when canonical code attempts to use a quarantined path."""


def validate_run_id(run_id: str) -> str:
    """Reject path separators and shell-significant text in a run identifier."""

    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", run_id):
        raise CanonicalPathError("Invalid run ID: use only letters, digits, dot, underscore, and hyphen")
    return run_id


def assert_canonical_path(path: Path | str, repository_root: Path | str) -> Path:
    """Resolve *path* and reject anything inside the noncanonical archive."""

    resolved = Path(path).resolve(strict=False)
    archive = (Path(repository_root) / "archive_noncanonical").resolve(strict=False)
    if resolved == archive or resolved.is_relative_to(archive):
        raise CanonicalPathError(f"Canonical analysis cannot use the noncanonical archive: {resolved}")
    return resolved


def repository_relative(path: Path | str, repository_root: Path | str) -> Path:
    """Normalize either a relative CLI path or an absolute repository path."""

    root = Path(repository_root).resolve(strict=False)
    candidate = Path(path)
    resolved = (root / candidate).resolve(strict=False) if not candidate.is_absolute() else candidate.resolve(strict=False)
    if resolved != root and not resolved.is_relative_to(root):
        raise CanonicalPathError(f"Path escapes repository root: {path}")
    assert_canonical_path(resolved, root)
    return resolved.relative_to(root)
