#!/usr/bin/env python3
"""Stream preprocessing and map every complete read pair to both organelles."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import tomllib
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

from organelle_pipeline.commands import MappingInputs, build_mapping_command
from organelle_pipeline.inventory import ACCEPTABLE_SOURCE_VALIDATION_STATUSES
from organelle_pipeline.paths import assert_canonical_path, repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    StaleOutputError,
    build_stage_fingerprint_from_hashes,
    mapping_pipeline_code_digest,
    sha256_file,
    sha256_json,
    validate_recorded_tool_versions,
    validate_resume,
    validate_saved_outputs,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--threads-per-job", type=int, default=8)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tool_versions() -> dict[str, str]:
    commands = {
        "fastp": ["fastp", "--version"],
        "bwa": ["bwa"],
        "samtools": ["samtools", "--version"],
    }
    result = {}
    for name, command in commands.items():
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        lines = [line.strip() for line in (completed.stdout + "\n" + completed.stderr).splitlines() if line.strip()]
        result[name] = " | ".join(lines[:3]) if lines else f"exit={completed.returncode}"
    return result


def pre_reconciliation_command(
    canonical_command: str,
    output_bam: Path,
    sample_id: str,
    maximum_unqualified_base_percent: int,
) -> str:
    """Reconstruct the exact earlier command for one-time provenance reconciliation."""

    command = canonical_command.replace(
        f"--unqualified_percent_limit {maximum_unqualified_base_percent} ",
        "",
    )
    sort_temp_dir = output_bam.parent / "sort_tmp"
    for label in ("name", "coordinate"):
        prefix = sort_temp_dir / f"{sample_id}.{label}"
        command = command.replace(f"-T {shlex.quote(str(prefix))} ", "")
    return command


def mapping_row_projection(row: dict[str, str]) -> dict[str, object]:
    return {
        "sample_id": row["sample_id"],
        "r1_paths": row["r1_paths"].split(";") if row["r1_paths"] else [],
        "r2_paths": row["r2_paths"].split(";") if row["r2_paths"] else [],
        "pair_status": row["pair_status"],
        "analysis_eligible": row["analysis_eligible"],
    }


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    if not config["mapping"]["mark_duplicates"]:
        raise RuntimeError("Canonical mapping requires duplicate marking/removal")
    if config["mapping"]["nuclear_decoy"]:
        raise RuntimeError("A nuclear decoy was configured but none is available in immutable sources")
    required_exclusions = 4 | 256 | 512 | 1024 | 2048
    if int(config["mapping"]["exclude_sam_flags"]) & required_exclusions != required_exclusions:
        raise RuntimeError("Canonical mapping must exclude unmapped, secondary, QC-failed, duplicate, and supplementary reads")
    manifest = assert_canonical_path(root / config["paths"]["sample_manifest"], root)
    source_manifest = root / "canonical_publication/provenance/manifests/source_inputs.tsv"
    reference = root / "canonical_publication/references/selected/organelle_combined.fa"
    assert_canonical_path(source_manifest, root)
    assert_canonical_path(reference, root)
    prerequisite_paths = {
        stage: root / "canonical_publication/provenance/runs" / args.run_id / f"{stage}.json"
        for stage in ("source_validation", "references", "metadata")
    }
    prerequisite_states = {}
    for stage, path in prerequisite_paths.items():
        if not path.is_file():
            raise RuntimeError(f"Mapping requires completed {stage} provenance")
        state = json.loads(path.read_text())
        acceptable_statuses = ACCEPTABLE_SOURCE_VALIDATION_STATUSES if stage == "source_validation" else {"complete"}
        if state.get("status") not in acceptable_statuses:
            raise RuntimeError(f"Mapping prerequisite {stage} did not complete successfully")
        validate_saved_outputs(root, state)
        prerequisite_states[stage] = state
    mapping_upstream = {stage: state["fingerprint"]["digest"] for stage, state in prerequisite_states.items()}
    source_hashes = {row["source_path"]: row["sha256"] for row in read_tsv(source_manifest)}
    samples = [row for row in read_tsv(manifest) if row["analysis_eligible"] == "yes" and row["pair_status"] == "complete"]
    if args.limit is not None:
        samples = samples[: args.limit]
    bam_dir = root / "canonical_publication/work" / args.run_id / "mapping"
    qc_dir = root / "canonical_publication/results/qc" / args.run_id / "fastp"
    log_dir = root / "canonical_publication/provenance/runs" / args.run_id / "logs/mapping"
    state_dir = root / "canonical_publication/provenance/runs" / args.run_id / "mapping"
    for directory in (bam_dir, qc_dir, log_dir, state_dir):
        directory.mkdir(parents=True, exist_ok=True)
    (bam_dir / "sort_tmp").mkdir(exist_ok=True)
    versions = tool_versions()
    git_commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    reference_hash = sha256_file(reference)
    config_hash = sha256_file(config_path)
    code_digest = mapping_pipeline_code_digest(root)
    current_executable_hashes = {}
    for tool in ("fastp", "bwa", "samtools"):
        executable = shutil.which(tool)
        if executable is None:
            raise RuntimeError(f"Current {tool} executable is unavailable")
        current_executable_hashes[tool] = sha256_file(Path(executable).resolve())

    def run_sample(row: dict[str, str]) -> str:
        sample_id = row["sample_id"]
        r1 = tuple(Path(value) for value in row["r1_paths"].split(";") if value)
        r2 = tuple(Path(value) for value in row["r2_paths"].split(";") if value)
        for path in (*r1, *r2):
            assert_canonical_path(root / path, root)
            if path.as_posix() not in source_hashes:
                raise RuntimeError(f"No immutable source checksum for {path}")
        bam = bam_dir / f"{sample_id}.organelle.bam"
        fastp_json = qc_dir / f"{sample_id}.fastp.json"
        fastp_html = qc_dir / f"{sample_id}.fastp.html"
        command = build_mapping_command(
            MappingInputs(
                sample_id=sample_id,
                r1_paths=r1,
                r2_paths=r2,
                reference=reference.relative_to(root),
                output_bam=bam.relative_to(root),
                fastp_json=fastp_json.relative_to(root),
                fastp_html=fastp_html.relative_to(root),
            ),
            threads=args.threads_per_job,
            qualified_quality_phred=int(config["preprocessing"]["qualified_quality_phred"]),
            maximum_unqualified_base_percent=int(config["preprocessing"]["maximum_unqualified_base_percent"]),
            minimum_length=int(config["preprocessing"]["minimum_length"]),
            minimum_mapping_quality=int(config["mapping"]["minimum_mapping_quality"]),
            exclude_sam_flags=int(config["mapping"]["exclude_sam_flags"]),
            detect_adapters_for_pe=bool(config["preprocessing"]["detect_adapters_for_pe"]),
        )
        declared_hashes = {
            **{path.as_posix(): source_hashes[path.as_posix()] for path in (*r1, *r2)},
            reference.relative_to(root).as_posix(): reference_hash,
            config_path.relative_to(root).as_posix(): config_hash,
        }
        fingerprint = build_stage_fingerprint_from_hashes(f"mapping:{sample_id}", declared_hashes, mapping_upstream, [command])
        state_path = state_dir / f"{sample_id}.json"
        if args.resume and state_path.exists():
            saved = json.loads(state_path.read_text())
            if saved.get("schema_version") == 2:
                rebind_reasons = []
                expected_policy = {
                    "preprocessing": config["preprocessing"],
                    "mapping": config["mapping"],
                    "threads_per_job": args.threads_per_job,
                }
                if saved.get("mapping_policy") != expected_policy or saved.get("current_reproduction_commands") != [command]:
                    raise StaleOutputError(f"Finalized mapping policy is stale for {sample_id}")
                finalized_inputs = dict(saved["fingerprint"]["inputs"])
                for input_path in (*r1, *r2):
                    if finalized_inputs.get(input_path.as_posix()) != source_hashes[input_path.as_posix()]:
                        raise StaleOutputError(f"Finalized read fingerprint is stale for {sample_id}: {input_path}")
                reference_key = reference.relative_to(root).as_posix()
                if finalized_inputs.get(reference_key) != reference_hash:
                    raise StaleOutputError(f"Finalized reference fingerprint is stale for {sample_id}")
                if finalized_inputs.get("metadata:mapping_row") != sha256_json(mapping_row_projection(row)):
                    raise StaleOutputError(f"Finalized metadata row is stale for {sample_id}")
                if finalized_inputs.get("provenance:pipeline_code_at_reconciliation") != code_digest:
                    rebind_reasons.append("loaded mapping code digest changed")
                validate_recorded_tool_versions(dict(saved.get("tool_versions", {})), versions)
                for tool, digest in current_executable_hashes.items():
                    if finalized_inputs.get(f"provenance:executable_sha256:{tool}") != digest:
                        raise StaleOutputError(f"Finalized {tool} executable is stale for {sample_id}")
                finalized_upstream = dict(saved["fingerprint"]["upstream"])
                for stage in ("source_validation", "references", "metadata"):
                    upstream_path = root / "canonical_publication/provenance/runs" / args.run_id / f"{stage}.json"
                    if not upstream_path.is_file():
                        raise StaleOutputError(f"Finalized mapping lacks current {stage} state for {sample_id}")
                    upstream_state = json.loads(upstream_path.read_text())
                    if finalized_upstream.get(stage) != upstream_state.get("fingerprint", {}).get("digest"):
                        rebind_reasons.append(f"{stage} fingerprint changed with identical mapping inputs")
            else:
                try:
                    validate_resume(saved["fingerprint"]["digest"], fingerprint)
                except StaleOutputError:
                    saved_inputs = dict(saved["fingerprint"]["inputs"])
                    config_key = config_path.relative_to(root).as_posix()
                    expected_nonconfig = {key: value for key, value in declared_hashes.items() if key != config_key}
                    observed_nonconfig = {key: value for key, value in saved_inputs.items() if key != config_key}
                    legacy_command = pre_reconciliation_command(
                        command,
                        bam.relative_to(root),
                        sample_id,
                        int(config["preprocessing"]["maximum_unqualified_base_percent"]),
                    )
                    recorded_fastp = str(saved.get("tool_versions", {}).get("fastp", ""))
                    if (
                        observed_nonconfig != expected_nonconfig
                        or saved["fingerprint"].get("commands") != [legacy_command]
                        or int(config["preprocessing"]["maximum_unqualified_base_percent"]) != 40
                        or not recorded_fastp.startswith("fastp 1.3.6")
                    ):
                        raise
            validate_saved_outputs(root, saved)
            suffix = (
                f"; pending provenance rebind: {', '.join(rebind_reasons)}" if saved.get("schema_version") == 2 and rebind_reasons else ""
            )
            return f"resume-valid {sample_id} (provisional or finalized provenance){suffix}"
        if bam.exists() or state_path.exists():
            raise RuntimeError(f"Existing incomplete or unvalidated mapping for {sample_id}; use a new run ID")
        log_path = log_dir / f"{sample_id}.log"
        with log_path.open("w") as log:
            completed = subprocess.run(
                ["bash", "-o", "pipefail", "-c", command],
                cwd=root,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                env=os.environ.copy(),
                check=False,
            )
        if completed.returncode != 0:
            raise RuntimeError(f"Mapping failed for {sample_id} with exit {completed.returncode}; see {log_path}")
        subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)
        outputs = {
            bam.relative_to(root).as_posix(): sha256_file(bam),
            bam.with_suffix(".bam.bai").relative_to(root).as_posix(): sha256_file(bam.with_suffix(".bam.bai")),
            fastp_json.relative_to(root).as_posix(): sha256_file(fastp_json),
            fastp_html.relative_to(root).as_posix(): sha256_file(fastp_html),
        }
        state = {
            "status": "complete",
            "sample_id": sample_id,
            "git_commit": git_commit,
            "tool_versions": versions,
            "fingerprint": asdict(fingerprint),
            "outputs": outputs,
        }
        state_path.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
        return f"mapped {sample_id}"

    failures = []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run_sample, row): row["sample_id"] for row in samples}
        for future in as_completed(futures):
            sample_id = futures[future]
            try:
                print(future.result(), flush=True)
            except Exception as error:
                failures.append((sample_id, str(error)))
                print(f"FAILED {sample_id}: {error}", flush=True)
    if failures:
        raise RuntimeError(f"{len(failures)} mappings failed: {failures[:3]}")
    print(f"completed or validated {len(samples)} sample mappings")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
