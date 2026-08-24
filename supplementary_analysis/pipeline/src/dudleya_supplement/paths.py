"""Filesystem boundaries for the supplementary workspace."""

from __future__ import annotations

from pathlib import Path


class SupplementPathError(ValueError):
    """Raised when a write would escape the supplementary workspace."""


def assert_output_path(path: Path | str, repository_root: Path | str) -> Path:
    root = Path(repository_root).resolve(strict=False)
    resolved = Path(path).resolve(strict=False)
    allowed = (root / "supplementary_analysis").resolve(strict=False)
    if resolved == allowed or resolved.is_relative_to(allowed):
        return resolved
    raise SupplementPathError(f"Supplementary output must remain beneath {allowed}: {resolved}")


def repository_input(path: Path | str, repository_root: Path | str) -> Path:
    root = Path(repository_root).resolve(strict=False)
    candidate = Path(path)
    resolved = (root / candidate).resolve(strict=False) if not candidate.is_absolute() else candidate.resolve(strict=False)
    if resolved != root and not resolved.is_relative_to(root):
        raise SupplementPathError(f"Input escapes repository root: {path}")
    archive = (root / "archive_noncanonical").resolve(strict=False)
    if resolved == archive or resolved.is_relative_to(archive):
        raise SupplementPathError(f"Archived input is forbidden: {resolved}")
    return resolved
