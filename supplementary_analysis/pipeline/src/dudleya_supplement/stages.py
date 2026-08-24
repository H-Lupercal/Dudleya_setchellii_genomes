"""Supplementary stage implementations."""

from __future__ import annotations

import json
import platform
import shutil
import subprocess
import tomllib
from dataclasses import asdict
from pathlib import Path

from .comparative_analysis import run_comparative_analyses
from .documentation import write_claim_decisions, write_inheritance_evidence, write_phase1_acceptance
from .finalization import write_acceptance, write_reports
from .identity_audit import run_identity_audit
from .io import read_tsv, write_json, write_tsv
from .likelihood import run_likelihood_mapping
from .metadata import apply_metadata_policy, derive_populations, verification_rows
from .provenance import (
    StaleSupplementError,
    build_fingerprint,
    code_input_hashes,
    filesystem_snapshot,
    sha256_file,
    validate_immutable_snapshot,
    validate_resume,
)
from .rendering import render_all_figures
from .scenario import run_all_sensitivity


def _git_commit(root: Path) -> str:
    return subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, check=True).stdout.strip()


def _state_path(root: Path, run_id: str, stage: str) -> Path:
    return root / "supplementary_analysis/provenance/runs" / run_id / f"{stage}.json"


def _load_state(root: Path, run_id: str, stage: str) -> dict[str, object]:
    path = _state_path(root, run_id, stage)
    if not path.is_file():
        raise RuntimeError(f"Missing upstream supplementary state: {path.relative_to(root)}")
    return json.loads(path.read_text())


def _finish_state(
    root: Path,
    run_id: str,
    stage: str,
    fingerprint: object,
    outputs: list[Path],
    extra: dict[str, object] | None = None,
) -> None:
    payload: dict[str, object] = {
        "status": "complete",
        "stage": stage,
        "fingerprint": asdict(fingerprint),
        "outputs": {path.relative_to(root).as_posix(): sha256_file(path) for path in outputs},
        "executable_versions": _software_versions(),
    }
    payload.update(extra or {})
    write_json(_state_path(root, run_id, stage), payload, root)


def _software_versions() -> dict[str, str]:
    versions = {"python": platform.python_version()}
    commands = {
        "git": ["git", "--version"],
        "bcftools": ["bcftools", "--version"],
        "samtools": ["samtools", "--version"],
        "iqtree3": ["iqtree3", "--version"],
        "mash": ["mash", "--version"],
        "outline": ["outline", "--help"],
    }
    for name, command in commands.items():
        if shutil.which(command[0]) is None:
            versions[name] = "unavailable"
            continue
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        versions[name] = (result.stdout or result.stderr).splitlines()[0].strip()
    return versions


def _canonical_fingerprint_value(state: dict[str, object]) -> str:
    """Normalize canonical state schemas used across the base run."""
    value = state.get("fingerprint")
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        digest = value.get("digest")
        return digest if isinstance(digest, str) else ""
    return ""


def _resume_or_fail(root: Path, run_id: str, stage: str, fingerprint: object, resume: bool) -> bool:
    path = _state_path(root, run_id, stage)
    if resume and path.is_file():
        saved = json.loads(path.read_text())
        validate_resume(saved["fingerprint"]["digest"], fingerprint)  # type: ignore[index,union-attr]
        for output, digest in saved["outputs"].items():  # type: ignore[union-attr]
            if sha256_file(root / output) != digest:
                raise StaleSupplementError(f"Changed supplementary output: {output}")
        print(f"resume-valid {stage}")
        return True
    if path.exists():
        raise RuntimeError(f"Existing supplementary state for {stage}; use --resume or a new run ID")
    return False


def verify_canonical_unchanged(root: Path, run_id: str) -> dict[str, dict[str, object]]:
    snapshot_path = root / "supplementary_analysis/provenance/manifests/publication-20260817.canonical_filesystem_snapshot.json"
    current = filesystem_snapshot(root / "canonical_publication")
    if snapshot_path.is_file():
        saved = json.loads(snapshot_path.read_text())
        validate_immutable_snapshot(saved, current)
    else:
        write_json(snapshot_path, current, root)
        blobs = subprocess.run(
            ["git", "ls-files", "-s", "canonical_publication"], cwd=root, capture_output=True, text=True, check=True
        ).stdout
        blob_path = root / "supplementary_analysis/provenance/manifests/publication-20260817.tracked_blobs.tsv"
        rows = []
        for line in blobs.splitlines():
            mode, digest, stage_path = line.split(maxsplit=2)
            stage_number, path = stage_path.split("\t", 1)
            rows.append({"path": path, "mode": mode, "git_blob": digest, "stage": stage_number})
        write_tsv(blob_path, rows, ["path", "mode", "git_blob", "stage"], root)
    fingerprint_path = root / "supplementary_analysis/provenance/manifests/publication-20260817.stage_fingerprints.tsv"
    if not fingerprint_path.is_file():
        fingerprint_rows = []
        canonical_run = root / "canonical_publication/provenance/runs/publication-20260817"
        for state_path in sorted(canonical_run.glob("*.json")):
            state = json.loads(state_path.read_text())
            digest = _canonical_fingerprint_value(state)
            fingerprint_rows.append(
                {
                    "path": state_path.relative_to(root).as_posix(),
                    "stage": state.get("stage", state_path.stem),
                    "fingerprint": digest,
                    "sha256": sha256_file(state_path),
                }
            )
        write_tsv(
            fingerprint_path,
            fingerprint_rows,
            ["path", "stage", "fingerprint", "sha256"],
            root,
        )
    acceptance = json.loads((root / "canonical_publication/provenance/runs/publication-20260817/ACCEPTANCE.json").read_text())
    if acceptance.get("status") != "PASS" or acceptance.get("sample_counts") != {
        "chloroplast": 276,
        "mitochondria": 271,
        "shared": 271,
    }:
        raise RuntimeError("Canonical base acceptance or sample counts do not match the approved baseline")
    return current


def run_canonical_guard(root: Path, config_path: Path, run_id: str, resume: bool, final: bool = False) -> None:
    config = tomllib.loads(config_path.read_text())
    canonical_config = config["canonical"]
    for path_key, hash_key in (
        ("acceptance", "acceptance_sha256"),
        ("sample_manifest", "sample_manifest_sha256"),
        ("chloroplast_alignment", "chloroplast_alignment_sha256"),
        ("mitochondria_alignment", "mitochondria_alignment_sha256"),
    ):
        path = root / canonical_config[path_key]
        observed = sha256_file(path)
        if observed != canonical_config[hash_key]:
            raise StaleSupplementError(f"Canonical expected hash mismatch: {canonical_config[path_key]}")
    current = verify_canonical_unchanged(root, run_id)
    stage = "canonical_guard_final" if final else "canonical_guard"
    baseline = root / "supplementary_analysis/provenance/manifests/publication-20260817.canonical_filesystem_snapshot.json"
    inputs = {
        **code_input_hashes(root),
        config_path.relative_to(root).as_posix(): sha256_file(config_path),
        baseline.relative_to(root).as_posix(): sha256_file(baseline),
        "canonical:filesystem_snapshot": __import__("hashlib").sha256(json.dumps(current, sort_keys=True).encode()).hexdigest(),
    }
    upstream = {}
    if final:
        acceptance_state = _load_state(root, run_id, "acceptance")
        upstream["acceptance"] = acceptance_state["fingerprint"]["digest"]  # type: ignore[index]
    fingerprint = build_fingerprint(
        stage, inputs, upstream, ["verify canonical filesystem is byte- and metadata-unchanged"], _git_commit(root)
    )
    if _resume_or_fail(root, run_id, stage, fingerprint, resume):
        return
    if final:
        diff = subprocess.run(
            ["git", "diff", "--exit-code", "main", "--", "canonical_publication"],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        )
        status = subprocess.run(
            ["git", "status", "--short", "--", "canonical_publication"],
            cwd=root,
            capture_output=True,
            text=True,
            check=True,
        )
        if diff.returncode or diff.stdout or diff.stderr or status.stdout or status.stderr:
            raise StaleSupplementError("Canonical Git diff/status is not clean against main")
    _finish_state(root, run_id, stage, fingerprint, [], {"canonical_entry_count": len(current)})


def run_metadata(root: Path, config_path: Path, run_id: str, resume: bool) -> None:
    verify_canonical_unchanged(root, run_id)
    canonical_state = _load_state(root, run_id, "canonical_guard")
    samples_path = root / "canonical_publication/metadata/samples/samples.tsv"
    populations_path = root / "canonical_publication/metadata/populations/populations.tsv"
    ambiguity_path = root / "canonical_publication/metadata/qc/publication-20260817/source_metadata_ambiguities.tsv"
    inputs = code_input_hashes(root)
    inputs.update(
        {path.relative_to(root).as_posix(): sha256_file(path) for path in (config_path, samples_path, populations_path, ambiguity_path)}
    )
    fingerprint = build_fingerprint(
        "metadata",
        inputs,
        {"canonical_guard": canonical_state["fingerprint"]["digest"]},  # type: ignore[index]
        ["apply decision-plan v2.5 sample-level metadata policy"],
        _git_commit(root),
    )
    if _resume_or_fail(root, run_id, "metadata", fingerprint, resume):
        return
    corrected, changes = apply_metadata_policy(read_tsv(samples_path))
    populations = derive_populations(corrected)
    sample_output = root / "supplementary_analysis/metadata/samples/samples.corrected-20260824.tsv"
    population_output = root / f"supplementary_analysis/metadata/populations/populations.{run_id}.tsv"
    correction_output = root / f"supplementary_analysis/metadata/qc/{run_id}/metadata_correction_manifest.tsv"
    verification_output = root / f"supplementary_analysis/metadata/qc/{run_id}/metadata_verification.tsv"
    write_tsv(sample_output, corrected, list(corrected[0]), root)
    write_tsv(population_output, populations, ["popcode", "species", "population_name"], root)
    write_tsv(
        correction_output,
        changes,
        [
            "sample_id",
            "old_popcode",
            "new_popcode_or_EXCLUDED",
            "evidence_source",
            "decision_author",
            "decision_date",
            "confidence_or_unresolved",
        ],
        root,
    )
    write_tsv(verification_output, verification_rows(), ["entity", "issue", "status", "action"], root)
    _finish_state(
        root,
        run_id,
        "metadata",
        fingerprint,
        [sample_output, population_output, correction_output, verification_output],
        {"sample_count": len(corrected), "population_inference_count": len(populations), "excluded_sample_count": len(changes)},
    )


CANONICAL_IMPORTS = (
    "canonical_publication/pipeline/src/organelle_pipeline/analysis.py",
    "canonical_publication/pipeline/src/organelle_pipeline/consensus.py",
    "canonical_publication/pipeline/src/organelle_pipeline/haplotypes.py",
    "canonical_publication/pipeline/src/organelle_pipeline/ordination.py",
    "canonical_publication/pipeline/src/organelle_pipeline/popgen.py",
    "canonical_publication/pipeline/src/organelle_pipeline/variants.py",
)


def _run_action_stage(
    root: Path,
    config_path: Path,
    run_id: str,
    stage: str,
    upstream_stage: str,
    resume: bool,
    commands: list[str],
    action,
) -> None:
    verify_canonical_unchanged(root, run_id)
    upstream_state = _load_state(root, run_id, upstream_stage)
    imported = [root / path for path in CANONICAL_IMPORTS]
    inputs = code_input_hashes(root, imported)
    inputs[config_path.relative_to(root).as_posix()] = sha256_file(config_path)
    baseline = root / "supplementary_analysis/provenance/manifests/publication-20260817.canonical_filesystem_snapshot.json"
    inputs[baseline.relative_to(root).as_posix()] = sha256_file(baseline)
    fingerprint = build_fingerprint(
        stage,
        inputs,
        {upstream_stage: upstream_state["fingerprint"]["digest"]},  # type: ignore[index]
        commands,
        _git_commit(root),
    )
    if _resume_or_fail(root, run_id, stage, fingerprint, resume):
        return
    outputs = list(action())
    missing = [path for path in outputs if not path.is_file()]
    if missing:
        raise RuntimeError(f"Stage {stage} declared missing outputs: {missing[:3]}")
    _finish_state(root, run_id, stage, fingerprint, sorted(set(outputs)))


def run_stage(stage: str, root: Path, config_path: Path, run_id: str, resume: bool) -> None:
    config = tomllib.loads(config_path.read_text())
    if stage == "canonical_guard":
        run_canonical_guard(root, config_path, run_id, resume)
    elif stage == "canonical_guard_final":
        run_canonical_guard(root, config_path, run_id, resume, final=True)
    elif stage == "metadata":
        run_metadata(root, config_path, run_id, resume)
    elif stage == "identity":
        _run_action_stage(
            root,
            config_path,
            run_id,
            stage,
            "metadata",
            resume,
            ["mash sketch -r -k 31 -s 100000", "bcftools query DP:AD"],
            lambda: run_identity_audit(root, run_id),
        )
    elif stage == "sensitivity":
        _run_action_stage(
            root,
            config_path,
            run_id,
            stage,
            "identity",
            resume,
            ["reuse canonical BAMs; no mapping", "bcftools haploid calls", "IQ-TREE fixed-seed trees", "9999-permutation Procrustes"],
            lambda: run_all_sensitivity(root, run_id, config),
        )
    elif stage == "claims":
        _run_action_stage(
            root,
            config_path,
            run_id,
            stage,
            "sensitivity",
            resume,
            ["map sensitivity status to claims"],
            lambda: write_claim_decisions(root, run_id),
        )
    elif stage == "inheritance":
        _run_action_stage(
            root,
            config_path,
            run_id,
            stage,
            "claims",
            resume,
            ["apply conservative inheritance wording gate"],
            lambda: write_inheritance_evidence(root, run_id),
        )
    elif stage == "phase1_gate":
        _run_action_stage(
            root,
            config_path,
            run_id,
            stage,
            "inheritance",
            resume,
            ["evaluate phase-1 blocking conditions"],
            lambda: write_phase1_acceptance(root, run_id),
        )
    elif stage == "likelihood_mapping":
        _run_action_stage(
            root,
            config_path,
            run_id,
            stage,
            "phase1_gate",
            resume,
            ["iqtree3 -lmap 100000 -n 0 with fixed model and seed", "conditional SplitsPy 0.0.10 NeighborNet"],
            lambda: run_likelihood_mapping(root, run_id, config),
        )
    elif stage == "comparative_analyses":
        _run_action_stage(
            root,
            config_path,
            run_id,
            stage,
            "likelihood_mapping",
            resume,
            [
                "support contraction and multifurcating unrooted RF",
                "1000 site draws",
                "1000 n=4 pi draws",
                "PC-QC permutation tests",
                "5 kb coordinate tracks",
            ],
            lambda: run_comparative_analyses(root, run_id),
        )
    elif stage == "figures":
        _run_action_stage(
            root,
            config_path,
            run_id,
            stage,
            "comparative_analyses",
            resume,
            ["render exactly six figure families as PNG PDF SVG", "render separate presentation replacements"],
            lambda: render_all_figures(root, run_id),
        )
    elif stage == "reports":
        _run_action_stage(
            root,
            config_path,
            run_id,
            stage,
            "figures",
            resume,
            ["write manuscript support and table manifests"],
            lambda: write_reports(root, run_id),
        )
    elif stage == "acceptance":
        _run_action_stage(
            root,
            config_path,
            run_id,
            stage,
            "reports",
            resume,
            ["checksum final artifacts and evaluate acceptance", "update CURRENT_RUN only on PASS"],
            lambda: write_acceptance(root, run_id, canonical_unchanged=True),
        )
    else:
        raise ValueError(f"Unknown supplementary stage: {stage}")
