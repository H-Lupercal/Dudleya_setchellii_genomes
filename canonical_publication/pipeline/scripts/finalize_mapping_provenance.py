#!/usr/bin/env python3
"""Validate every mapping artifact and finalize dependency-complete provenance."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import tomllib
from dataclasses import asdict
from pathlib import Path

from organelle_pipeline.commands import MappingInputs, build_mapping_command
from organelle_pipeline.inventory import ACCEPTABLE_SOURCE_VALIDATION_STATUSES
from organelle_pipeline.paths import assert_canonical_path, repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    mapping_pipeline_code_digest,
    sha256_file,
    sha256_json,
    validate_recorded_tool_versions,
    validate_saved_outputs,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--threads-per-job", type=int, default=8)
    parser.add_argument("--resume", action="store_true")
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
    manifest = assert_canonical_path(root / config["paths"]["sample_manifest"], root)
    source_manifest = root / "canonical_publication/provenance/manifests/source_inputs.tsv"
    reference = root / "canonical_publication/references/selected/organelle_combined.fa"
    run_state_dir = root / "canonical_publication/provenance/runs" / args.run_id
    state_dir = run_state_dir / "mapping"
    completion_path = run_state_dir / "mapping_complete.json"
    bam_dir = root / "canonical_publication/work" / args.run_id / "mapping"
    fastp_dir = root / "canonical_publication/results/qc" / args.run_id / "fastp"
    prerequisite_paths = {
        "source_validation": run_state_dir / "source_validation.json",
        "references": run_state_dir / "references.json",
        "metadata": run_state_dir / "metadata.json",
    }
    for label, path in prerequisite_paths.items():
        if not path.is_file():
            raise RuntimeError(f"Mapping provenance requires completed {label} state: {path}")
    prerequisite_states = {label: json.loads(path.read_text()) for label, path in prerequisite_paths.items()}
    for label, state in prerequisite_states.items():
        acceptable_statuses = ACCEPTABLE_SOURCE_VALIDATION_STATUSES if label == "source_validation" else {"complete"}
        if state.get("status") not in acceptable_statuses:
            raise RuntimeError(f"Mapping provenance prerequisite did not pass: {label}")

    samples = [row for row in read_tsv(manifest) if row["analysis_eligible"] == "yes" and row["pair_status"] == "complete"]
    source_hashes = {row["source_path"]: row["sha256"] for row in read_tsv(source_manifest)}
    source_validation_output = root / next(iter(prerequisite_states["source_validation"]["outputs"]))
    provider_status_by_source = {
        row["resolved_source_path"]: row["status"] for row in read_tsv(source_validation_output) if row["resolved_source_path"]
    }
    current_versions = tool_versions()
    executable_hashes = {}
    for tool in ("fastp", "bwa", "samtools"):
        executable = shutil.which(tool)
        if executable is None:
            raise RuntimeError(f"Required mapping executable is unavailable: {tool}")
        executable_hashes[tool] = sha256_file(Path(executable).resolve())
    reference_hash = sha256_file(reference)
    mapping_policy = {
        "preprocessing": config["preprocessing"],
        "mapping": config["mapping"],
        "threads_per_job": args.threads_per_job,
    }
    policy_hash = sha256_json(mapping_policy)
    reconciliation_code_digest = mapping_pipeline_code_digest(root)
    reconciliation_git_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=root,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    upstream = {label: state["fingerprint"]["digest"] for label, state in prerequisite_states.items()}
    finalized: list[tuple[Path, dict[str, object]]] = []
    sample_fingerprints: dict[str, str] = {}
    for row in samples:
        sample_id = row["sample_id"]
        state_path = state_dir / f"{sample_id}.json"
        if not state_path.is_file():
            raise RuntimeError(f"Missing provisional mapping state for {sample_id}")
        saved = json.loads(state_path.read_text())
        if saved.get("status") != "complete" or saved.get("sample_id") != sample_id:
            raise RuntimeError(f"Invalid provisional mapping state for {sample_id}")
        validate_saved_outputs(root, saved)
        old_fingerprint = saved.get("pre_reconciliation_fingerprint", saved["fingerprint"])
        old_inputs = dict(old_fingerprint["inputs"])
        r1 = tuple(Path(value) for value in row["r1_paths"].split(";") if value)
        r2 = tuple(Path(value) for value in row["r2_paths"].split(";") if value)
        for read_path in (*r1, *r2):
            assert_canonical_path(root / read_path, root)
            expected_hash = source_hashes.get(read_path.as_posix())
            if expected_hash is None or old_inputs.get(read_path.as_posix()) != expected_hash:
                raise RuntimeError(f"Provisional input checksum mismatch for {sample_id}: {read_path}")
            if provider_status_by_source.get(read_path.as_posix()) != "PASS":
                raise RuntimeError(f"Mapped read lacks a passing provider checksum for {sample_id}: {read_path}")
        reference_key = reference.relative_to(root).as_posix()
        if old_inputs.get(reference_key) != reference_hash:
            raise RuntimeError(f"Provisional reference checksum mismatch for {sample_id}")

        bam = bam_dir / f"{sample_id}.organelle.bam"
        bai = bam.with_suffix(".bam.bai")
        fastp_json = fastp_dir / f"{sample_id}.fastp.json"
        fastp_html = fastp_dir / f"{sample_id}.fastp.html"
        for output in (bam, bai, fastp_json, fastp_html):
            if not output.is_file():
                raise RuntimeError(f"Missing mapping output for {sample_id}: {output}")
        subprocess.run(["samtools", "quickcheck", "-v", str(bam)], check=True)
        idxstats_lines = subprocess.run(
            ["samtools", "idxstats", str(bam)],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.splitlines()
        indexed_references = {line.split("\t", 1)[0] for line in idxstats_lines if line and not line.startswith("*\t")}
        if indexed_references != {"chloroplast", "mitochondria"}:
            raise RuntimeError(f"BAM index/reference validation failed for {sample_id}: {sorted(indexed_references)}")
        total_records = int(
            subprocess.run(
                ["samtools", "view", "-c", str(bam)],
                capture_output=True,
                text=True,
                check=True,
            ).stdout
        )
        forbidden_records = int(
            subprocess.run(
                ["samtools", "view", "-c", "--incl-flags", str(config["mapping"]["exclude_sam_flags"]), str(bam)],
                capture_output=True,
                text=True,
                check=True,
            ).stdout
        )
        passing_mapq_records = int(
            subprocess.run(
                ["samtools", "view", "-c", "-q", str(config["mapping"]["minimum_mapping_quality"]), str(bam)],
                capture_output=True,
                text=True,
                check=True,
            ).stdout
        )
        if forbidden_records or passing_mapq_records != total_records:
            raise RuntimeError(
                f"BAM filter validation failed for {sample_id}: total={total_records}, "
                f"forbidden={forbidden_records}, MAPQ-passing={passing_mapq_records}"
            )
        header = subprocess.run(
            ["samtools", "view", "-H", str(bam)],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        header_line = next((line for line in header.splitlines() if line.startswith("@HD\t")), "")
        header_fields = dict(field.split(":", 1) for field in header_line.split("\t")[1:] if ":" in field)
        if header_fields.get("SO") != "coordinate":
            raise RuntimeError(f"BAM coordinate-sort validation failed for {sample_id}: {header_fields}")
        read_groups = []
        for line in header.splitlines():
            if not line.startswith("@RG\t"):
                continue
            read_groups.append(dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field))
        if read_groups != [{"ID": sample_id, "SM": sample_id, "PL": "ILLUMINA"}]:
            raise RuntimeError(f"BAM read-group validation failed for {sample_id}: {read_groups}")

        canonical_command = build_mapping_command(
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
        executed_commands = saved.get("executed_commands", old_fingerprint["commands"])
        if len(executed_commands) != 1:
            raise RuntimeError(f"Expected one executed mapping command for {sample_id}")
        executed = executed_commands[0]
        required_fragments = (
            f"--qualified_quality_phred {config['preprocessing']['qualified_quality_phred']}",
            f"--length_required {config['preprocessing']['minimum_length']}",
            f"samtools view -F {config['mapping']['exclude_sam_flags']} -q {config['mapping']['minimum_mapping_quality']}",
            "samtools fixmate -m",
            "samtools markdup -r",
            f"ID:{sample_id}\\tSM:{sample_id}",
        )
        if any(fragment not in executed for fragment in required_fragments):
            raise RuntimeError(f"Executed mapping policy mismatch for {sample_id}")
        unqualified_fragment = f"--unqualified_percent_limit {config['preprocessing']['maximum_unqualified_base_percent']}"
        resolved_defaults = []
        if unqualified_fragment not in executed:
            if int(config["preprocessing"]["maximum_unqualified_base_percent"]) != 40 or not current_versions["fastp"].startswith(
                "fastp 1.3.6"
            ):
                raise RuntimeError(f"Cannot resolve fastp default policy for {sample_id}")
            resolved_defaults.append("fastp --unqualified_percent_limit=40 (pinned fastp 1.3.6 default)")

        mapping_commit = str(saved.get("git_commit") or "")
        if not mapping_commit:
            raise RuntimeError(f"Mapping Git commit was not captured for {sample_id}")
        mapping_versions, version_reconciliation = validate_recorded_tool_versions(
            dict(saved.get("tool_versions", {})),
            current_versions,
            # States produced by the already-running pre-reconciliation mapper
            # missed BWA's stderr-only banner.  No other missing version is
            # accepted, and the resolved executable hash is fingerprinted.
            reconcilable_missing={"bwa"},
        )
        resolved_defaults.extend(version_reconciliation)
        declared = {
            **{read_path.as_posix(): source_hashes[read_path.as_posix()] for read_path in (*r1, *r2)},
            reference_key: reference_hash,
            "configuration:mapping_policy": policy_hash,
            "metadata:mapping_row": sha256_json(mapping_row_projection(row)),
            "provenance:git_commit": mapping_commit,
            "provenance:reconciliation_git_commit": reconciliation_git_commit,
            "provenance:pipeline_code_at_reconciliation": reconciliation_code_digest,
            **{f"provenance:tool:{tool}": version for tool, version in mapping_versions.items()},
            **{f"provenance:executable_sha256:{tool}": digest for tool, digest in executable_hashes.items()},
        }
        fingerprint = build_stage_fingerprint_from_hashes(
            f"mapping:{sample_id}",
            declared,
            upstream,
            list(executed_commands),
        )
        outputs = {path.relative_to(root).as_posix(): sha256_file(path) for path in (bam, bai, fastp_json, fastp_html)}
        finalized_state: dict[str, object] = {
            "schema_version": 2,
            "status": "complete",
            "sample_id": sample_id,
            "git_commit": mapping_commit,
            "reconciliation_git_commit": reconciliation_git_commit,
            "tool_versions": mapping_versions,
            "mapping_policy": mapping_policy,
            "executed_commands": executed_commands,
            "current_reproduction_commands": [canonical_command],
            "resolved_defaults": resolved_defaults,
            "bam_validation": {
                "record_count": total_records,
                "forbidden_flag_record_count": forbidden_records,
                "mapq_passing_record_count": passing_mapq_records,
                "coordinate_sort_order": header_fields["SO"],
                "indexed_references": sorted(indexed_references),
                "read_group": read_groups[0],
            },
            "pre_reconciliation_fingerprint": old_fingerprint,
            "fingerprint": asdict(fingerprint),
            "outputs": outputs,
        }
        if saved.get("schema_version") == 2:
            validate_saved_outputs(root, saved)
            previous_digest = saved["fingerprint"]["digest"]
            if previous_digest == fingerprint.digest:
                finalized_state = saved
            else:
                finalized_state["provenance_rebind"] = {
                    "rebound_from_fingerprint": previous_digest,
                    "reason": (
                        "stage-scoped provenance refresh after all raw-read, reference, metadata, policy, "
                        "tool, executable, command, BAM, BAI, and fastp validations passed"
                    ),
                }
        finalized.append((state_path, finalized_state))
        sample_fingerprints[sample_id] = fingerprint.digest

    completion_fingerprint = build_stage_fingerprint_from_hashes(
        "mapping_provenance_completion",
        {
            "configuration:mapping_policy": policy_hash,
            "metadata:sample_manifest": sha256_file(manifest),
            "reference:combined": reference_hash,
        },
        {**upstream, **sample_fingerprints},
        ["validate all mapping inputs, exact executed commands, BAM/BAI/fastp outputs, and resolved pinned defaults"],
    )
    completion_rebind: dict[str, str] | None = None
    if args.resume and completion_path.exists():
        saved_completion = json.loads(completion_path.read_text())
        validate_saved_outputs(root, saved_completion)
        previous_digest = saved_completion["fingerprint"]["digest"]
        if previous_digest == completion_fingerprint.digest:
            print(f"resume-valid mapping provenance for {len(samples)} samples")
            return 0
        completion_rebind = {
            "rebound_from_fingerprint": previous_digest,
            "reason": "validated stage-scoped provenance refresh; no mapping artifact was regenerated",
        }
    if completion_path.exists() and completion_rebind is None:
        raise RuntimeError("Mapping completion state already exists; use --resume or a new run ID")
    for state_path, state in finalized:
        state_path.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
    completion_outputs = {state_path.relative_to(root).as_posix(): sha256_file(state_path) for state_path, _ in finalized}
    completion_state: dict[str, object] = {
        "status": "complete",
        "sample_count": len(samples),
        "tool_versions": current_versions,
        "reconciliation_git_commit": reconciliation_git_commit,
        "executable_sha256": executable_hashes,
        "fingerprint": asdict(completion_fingerprint),
        "outputs": completion_outputs,
    }
    if completion_rebind is not None:
        completion_state["provenance_rebind"] = completion_rebind
    completion_path.write_text(json.dumps(completion_state, indent=2, sort_keys=True) + "\n")
    print(f"finalized dependency-complete mapping provenance for {len(samples)} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
