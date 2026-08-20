#!/usr/bin/env python3
"""Verify immutable FASTQs against sequencing-provider MD5 manifests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import asdict
from pathlib import Path

from organelle_pipeline.inventory import classify_provider_md5, source_validation_status
from organelle_pipeline.paths import validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def file_digests(path: Path) -> tuple[str, str]:
    md5_digest = hashlib.md5(usedforsecurity=False)
    sha256_digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            md5_digest.update(block)
            sha256_digest.update(block)
    return md5_digest.hexdigest(), sha256_digest.hexdigest()


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    raw_root = root / "source_data/raw_reads"
    source_manifest = root / "canonical_publication/provenance/manifests/source_inputs.tsv"
    provider_manifests = sorted(raw_root.rglob("md5sum.txt"))
    if not provider_manifests:
        raise RuntimeError("No sequencing-provider MD5 manifests were found in immutable sources")
    with source_manifest.open(newline="") as handle:
        source_rows = list(csv.DictReader(handle, delimiter="\t"))
    source_by_path: dict[str, dict[str, str]] = {}
    for row in source_rows:
        source_path = row["source_path"]
        if source_path in source_by_path:
            raise RuntimeError(f"Duplicate immutable-source manifest path: {source_path}")
        source_by_path[source_path] = row
    observed_source_paths = {
        path.relative_to(root).as_posix() for path in (root / "source_data").rglob("*") if path.is_file() or path.is_symlink()
    }
    if observed_source_paths != set(source_by_path):
        missing = sorted(set(source_by_path) - observed_source_paths)
        unmanifested = sorted(observed_source_paths - set(source_by_path))
        raise RuntimeError(f"Immutable-source manifest/filesystem mismatch; missing={missing[:5]}, unmanifested={unmanifested[:5]}")
    output = root / "canonical_publication/provenance/manifests" / f"{args.run_id}.provider_md5_validation.tsv"
    state_path = root / "canonical_publication/provenance/runs" / args.run_id / "source_validation.json"
    state_path.parent.mkdir(parents=True, exist_ok=True)
    fingerprint = build_stage_fingerprint_from_hashes(
        "immutable_source_validation",
        {
            **runtime_provenance(root, {"python": ("python", "--version")}),
            source_manifest.relative_to(root).as_posix(): sha256_file(source_manifest),
            **{path.relative_to(root).as_posix(): sha256_file(path) for path in provider_manifests},
        },
        {},
        ["resolve each provider MD5 basename uniquely beneath its deposit and compare the FASTQ content digest"],
    )
    if args.resume and state_path.exists():
        saved = json.loads(state_path.read_text())
        validate_resume(saved["fingerprint"]["digest"], fingerprint)
        validate_saved_outputs(root, saved)
        print("resume fingerprint valid; re-reading immutable sources to confirm current MD5 content")
    elif state_path.exists() or output.exists():
        raise RuntimeError("Existing unvalidated source-MD5 output")
    rows = []
    failures = []
    missing_entries = []
    self_reference_warnings = []
    observed_sha256: dict[str, str] = {}
    source_inventory_passes = 0
    for provider_manifest in provider_manifests:
        deposit_root = provider_manifest.parent
        for line_number, line in enumerate(provider_manifest.read_text().splitlines(), 1):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 2:
                failures.append(f"malformed {provider_manifest}:{line_number}")
                continue
            expected, provider_name = fields
            provider_name = provider_name.lstrip("*")
            matches = sorted(deposit_root.rglob(Path(provider_name).name))
            if len(matches) == 0:
                missing_entries.append(provider_name)
                rows.append(
                    {
                        "provider_manifest": provider_manifest.relative_to(root).as_posix(),
                        "provider_name": provider_name,
                        "resolved_source_path": "",
                        "expected_md5": expected,
                        "observed_md5": "",
                        "status": "DECLARED_MISSING",
                    }
                )
                print(f"provider-md5 DECLARED_MISSING {provider_name}", flush=True)
                continue
            if len(matches) != 1:
                failures.append(f"{provider_manifest}:{line_number} resolved {len(matches)} files")
                rows.append(
                    {
                        "provider_manifest": provider_manifest.relative_to(root).as_posix(),
                        "provider_name": provider_name,
                        "resolved_source_path": "",
                        "expected_md5": expected,
                        "observed_md5": "",
                        "status": "FAIL_RESOLUTION",
                    }
                )
                continue
            source = matches[0]
            observed, source_sha256 = file_digests(source)
            source_key = source.relative_to(root).as_posix()
            observed_sha256[source_key] = source_sha256
            source_record = source_by_path[source_key]
            source_inventory_matches = (
                source_record["artifact_type"] == "file"
                and int(source_record["size_bytes"]) == source.stat().st_size
                and source_record["sha256"] == source_sha256
            )
            source_inventory_passes += int(source_inventory_matches)
            status = classify_provider_md5(
                expected=expected,
                observed=observed,
                source_is_provider_manifest=source.resolve() == provider_manifest.resolve(),
                source_inventory_matches=source_inventory_matches,
            )
            if status == "FAIL_CHECKSUM":
                failures.append(source.relative_to(root).as_posix())
            elif status == "UNVERIFIABLE_SELF_REFERENCE":
                self_reference_warnings.append(source.relative_to(root).as_posix())
            rows.append(
                {
                    "provider_manifest": provider_manifest.relative_to(root).as_posix(),
                    "provider_name": provider_name,
                    "resolved_source_path": source.relative_to(root).as_posix(),
                    "expected_md5": expected.lower(),
                    "observed_md5": observed.lower(),
                    "status": status,
                }
            )
            print(f"provider-md5 {status} {source.name}", flush=True)
    for source_key, source_record in source_by_path.items():
        if source_key in observed_sha256:
            continue
        source = root / source_key
        if source_record["artifact_type"] != "file" or source.is_symlink():
            failures.append(f"unsupported immutable-source artifact type: {source_key}")
            continue
        source_sha256 = sha256_file(source)
        observed_sha256[source_key] = source_sha256
        source_inventory_matches = int(source_record["size_bytes"]) == source.stat().st_size and source_record["sha256"] == source_sha256
        source_inventory_passes += int(source_inventory_matches)
        if not source_inventory_matches:
            failures.append(f"immutable-source SHA-256 mismatch: {source_key}")
    with output.open("w", newline="") as handle:
        fields = [
            "provider_manifest",
            "provider_name",
            "resolved_source_path",
            "expected_md5",
            "observed_md5",
            "status",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    state_path.write_text(
        json.dumps(
            {
                "status": source_validation_status(
                    has_failures=bool(failures),
                    has_declared_missing=bool(missing_entries),
                    has_self_reference_warning=bool(self_reference_warnings),
                ),
                "provider_manifest_entries": len(rows),
                "checksum_validated_files": sum(row["status"] == "PASS" for row in rows),
                "source_inventory_files": len(source_by_path),
                "source_inventory_sha256_pass": source_inventory_passes,
                "declared_missing_provider_entries": missing_entries,
                "unverifiable_provider_manifest_self_references": self_reference_warnings,
                "failures": failures,
                "fingerprint": asdict(fingerprint),
                "outputs": {output.relative_to(root).as_posix(): sha256_file(output)},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    if failures:
        raise RuntimeError(f"Source MD5 validation failed: {failures[:5]}")
    print(f"validated {sum(row['status'] == 'PASS' for row in rows)} provider checksums; declared missing entries={len(missing_entries)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
