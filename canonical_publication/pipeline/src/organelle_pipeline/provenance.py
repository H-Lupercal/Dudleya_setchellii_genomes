"""Content fingerprints and strict resume validation."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


class StaleOutputError(RuntimeError):
    """Raised when a saved stage fingerprint no longer matches its inputs."""


@dataclass(frozen=True)
class StageFingerprint:
    stage: str
    inputs: tuple[tuple[str, str], ...]
    upstream: tuple[tuple[str, str], ...]
    commands: tuple[str, ...]
    digest: str


def sha256_file(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_json(value: object) -> str:
    """Hash a JSON-compatible value using a deterministic serialization."""

    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _code_digest(repository_root: Path | str, files: list[Path]) -> str:
    root = Path(repository_root).resolve()
    digest = hashlib.sha256()
    for path in sorted(files):
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(b"\0")
        digest.update(sha256_file(path).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def pipeline_code_digest(repository_root: Path | str) -> str:
    """Hash every canonical implementation file for repository-wide audits."""

    root = Path(repository_root).resolve()
    pipeline = root / "canonical_publication/pipeline"
    files = [
        path
        for path in pipeline.rglob("*")
        if path.is_file() and path.suffix in {".py", ".sh", ".toml"} and "__pycache__" not in path.parts
    ]
    return _code_digest(root, files)


def runtime_pipeline_code_digest(
    repository_root: Path | str,
    loaded_files: tuple[Path, ...] | list[Path] | None = None,
) -> str:
    """Hash only code loaded by the executing stage, preserving DAG direction."""

    root = Path(repository_root).resolve()
    pipeline = (root / "canonical_publication/pipeline").resolve()
    candidates = loaded_files
    if candidates is None:
        candidates = tuple(
            Path(module_file) for module in sys.modules.values() if isinstance((module_file := getattr(module, "__file__", None)), str)
        )
    files = []
    for candidate in candidates:
        path = Path(candidate).resolve()
        if path.is_file() and path.suffix == ".py" and path.is_relative_to(pipeline):
            files.append(path)
    if not files:
        # Tests and embedded callers may not execute from a canonical script,
        # but the provenance implementation itself is always stage-relevant.
        provenance_file = Path(__file__).resolve()
        if provenance_file.is_relative_to(pipeline):
            files.append(provenance_file)
    return _code_digest(root, sorted(set(files)))


def mapping_pipeline_code_digest(repository_root: Path | str) -> str:
    """Hash the mapping executor's fixed code surface for later validation.

    Mapping provenance is finalized by a separate entry script.  A runtime
    module scan there would hash the validator rather than the executor and
    make every successful resume immediately stale again.  Both entrypoints
    therefore use this same explicit mapping-stage code surface.
    """

    root = Path(repository_root).resolve()
    pipeline = root / "canonical_publication/pipeline"
    relative_files = (
        "scripts/map_samples.py",
        "src/organelle_pipeline/__init__.py",
        "src/organelle_pipeline/commands.py",
        "src/organelle_pipeline/inventory.py",
        "src/organelle_pipeline/metadata.py",
        "src/organelle_pipeline/paths.py",
        "src/organelle_pipeline/provenance.py",
    )
    files = [pipeline / relative for relative in relative_files]
    missing = [path for path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Mapping code surface is incomplete: {missing}")
    return _code_digest(root, files)


def _command_version(command: tuple[str, ...]) -> str:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    lines = [line.strip() for line in (completed.stdout + "\n" + completed.stderr).splitlines() if line.strip()]
    if not lines:
        return f"exit={completed.returncode}; no version text"
    return " | ".join(lines[:3])


def runtime_provenance(
    repository_root: Path | str,
    tool_commands: dict[str, tuple[str, ...]],
) -> dict[str, str]:
    """Return readable execution facts that participate in a stage digest."""

    root = Path(repository_root).resolve()
    git_commit = _command_version(("git", "rev-parse", "HEAD"))
    recorded = {
        "provenance:git_commit": git_commit,
        "provenance:pipeline_code": runtime_pipeline_code_digest(root),
    }
    for name, command in sorted(tool_commands.items()):
        recorded[f"provenance:tool:{name}"] = _command_version(command)
        executable = shutil.which(command[0])
        recorded[f"provenance:executable_sha256:{name}"] = sha256_file(Path(executable).resolve()) if executable else "unresolved"
    return recorded


def build_stage_fingerprint(
    stage: str,
    inputs: list[Path] | tuple[Path, ...],
    upstream: dict[str, str],
    commands: list[str] | tuple[str, ...],
) -> StageFingerprint:
    input_hashes = {str(Path(path)): sha256_file(path) for path in inputs}
    return build_stage_fingerprint_from_hashes(stage, input_hashes, upstream, commands)


def build_stage_fingerprint_from_hashes(
    stage: str,
    input_hashes: dict[str, str],
    upstream: dict[str, str],
    commands: list[str] | tuple[str, ...],
) -> StageFingerprint:
    """Fingerprint a stage using hashes from a verified immutable manifest."""

    normalized_inputs = tuple(sorted(input_hashes.items()))
    upstream_hashes = tuple(sorted(upstream.items()))
    command_tuple = tuple(commands)
    payload = {
        "stage": stage,
        "inputs": normalized_inputs,
        "upstream": upstream_hashes,
        "commands": command_tuple,
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return StageFingerprint(digest=digest, **payload)


def validate_resume(saved_digest: str, current: StageFingerprint) -> None:
    if saved_digest != current.digest:
        raise StaleOutputError(f"Stage {current.stage} is stale: saved {saved_digest}, current {current.digest}")


def validate_saved_outputs(repository_root: Path | str, saved_state: dict[str, object]) -> None:
    """Reject resume when any declared output is missing or content-changed."""

    root = Path(repository_root).resolve()
    outputs = saved_state.get("outputs")
    if not isinstance(outputs, dict):
        raise StaleOutputError("Saved state has no output checksum mapping")
    for relative, expected in outputs.items():
        path = (root / str(relative)).resolve()
        try:
            path.relative_to(root)
        except ValueError as error:
            raise StaleOutputError(f"Saved output escapes repository: {relative}") from error
        if not path.is_file():
            raise StaleOutputError(f"Saved output is missing: {relative}")
        observed = sha256_file(path)
        if observed != expected:
            raise StaleOutputError(f"Saved output checksum mismatch: {relative}")


def validate_recorded_tool_versions(
    recorded: dict[str, str],
    current: dict[str, str],
    reconcilable_missing: set[str] | frozenset[str] = frozenset(),
) -> tuple[dict[str, str], tuple[str, ...]]:
    """Validate execution-tool versions, explicitly noting legacy reconciliation.

    Earlier provisional mapping states did not capture BWA's stderr-only version
    banner.  A missing version may be reconciled only when the caller names that
    tool explicitly; the pinned current executable and its SHA-256 must then be
    included in the finalized fingerprint.
    """

    resolved: dict[str, str] = {}
    notes: list[str] = []
    for tool, current_version in sorted(current.items()):
        if not current_version:
            raise StaleOutputError(f"Current {tool} version could not be resolved")
        recorded_version = recorded.get(tool, "")
        if not recorded_version:
            if tool not in reconcilable_missing:
                raise StaleOutputError(f"Recorded {tool} version is missing")
            notes.append(f"{tool} version absent in provisional state; resolved from pinned executable during reconciliation")
            resolved[tool] = current_version
            continue
        if not current_version.startswith(recorded_version):
            raise StaleOutputError(
                f"Current {tool} version differs from recorded execution version: {current_version!r} != {recorded_version!r}"
            )
        resolved[tool] = current_version
    return resolved, tuple(notes)


def write_fingerprint(path: Path | str, fingerprint: StageFingerprint) -> None:
    Path(path).write_text(json.dumps(asdict(fingerprint), indent=2, sort_keys=True) + "\n")
