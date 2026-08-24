"""Content-addressed provenance for supplementary stages."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from pathlib import Path


class StaleSupplementError(RuntimeError):
    """Raised when a saved stage cannot be reused safely."""


@dataclass(frozen=True)
class StageFingerprint:
    stage: str
    digest: str
    inputs: dict[str, str]
    upstream: dict[str, str]
    commands: tuple[str, ...]
    git_commit: str


def sha256_file(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_fingerprint(
    stage: str,
    inputs: dict[str, str],
    upstream: dict[str, str],
    commands: list[str] | tuple[str, ...],
    git_commit: str,
) -> StageFingerprint:
    payload = {
        "stage": stage,
        "inputs": dict(sorted(inputs.items())),
        "upstream": dict(sorted(upstream.items())),
        "commands": list(commands),
        "git_commit": git_commit,
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return StageFingerprint(stage, digest, payload["inputs"], payload["upstream"], tuple(commands), git_commit)  # type: ignore[arg-type]


def validate_resume(saved_digest: str, current: StageFingerprint) -> None:
    if saved_digest != current.digest:
        raise StaleSupplementError(f"Supplementary stage {current.stage} is stale: saved {saved_digest}, current {current.digest}")


def filesystem_snapshot(root: Path | str) -> dict[str, dict[str, object]]:
    base = Path(root)
    result: dict[str, dict[str, object]] = {}
    for path in sorted(base.rglob("*")):
        stat = path.lstat()
        relative = path.relative_to(base).as_posix()
        if path.is_symlink():
            kind = "symlink"
            target: str | None = os.readlink(path)
        elif path.is_dir():
            kind = "directory"
            target = None
        else:
            kind = "file"
            target = None
        result[relative] = {
            "type": kind,
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "target": target,
            "sha256": sha256_file(path) if kind == "file" else None,
        }
    return result


def validate_immutable_snapshot(saved: dict[str, dict[str, object]], current: dict[str, dict[str, object]]) -> None:
    if saved != current:
        changed = sorted(set(saved) ^ set(current))
        if not changed:
            changed = [key for key in saved if saved[key] != current[key]]
        preview = ", ".join(changed[:5])
        raise StaleSupplementError(f"immutable canonical filesystem changed: {preview}")


def code_input_hashes(repository_root: Path | str, imported_canonical_files: list[Path] | tuple[Path, ...] = ()) -> dict[str, str]:
    root = Path(repository_root).resolve()
    files = sorted((root / "supplementary_analysis/pipeline/src").rglob("*.py"))
    files.extend(sorted((root / "supplementary_analysis/pipeline/scripts").glob("*.py")))
    files.extend(Path(path).resolve() for path in imported_canonical_files)
    result: dict[str, str] = {}
    for path in files:
        relative = path.relative_to(root).as_posix()
        result[relative] = sha256_file(path)
    return dict(sorted(result.items()))
