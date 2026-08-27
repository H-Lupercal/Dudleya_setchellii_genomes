import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from organelle_pipeline.provenance import sha256_file
from organelle_pipeline.publication_package import (
    PublicationPackageError,
    create_publication_package,
    validate_package_manifest,
    verify_publication_package,
)

BASE_RUN_ID = "test-run"
PACKAGE_ID = "test-run-package"


def bytes_digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def write_tsv(path: Path, rows: list[dict[str, str]], fields: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def miniature_publication(root: Path) -> None:
    unchanged_paths = {
        "canonical_publication/results/run/result.tsv": b"result\n1\n",
        "canonical_publication/metadata/samples/samples.tsv": b"sample\nDU1\n",
        "canonical_publication/references/selected/reference.fa": b">ref\nACGT\n",
        f"canonical_publication/provenance/runs/{BASE_RUN_ID}/qc.json": b'{"status":"complete"}\n',
        "canonical_publication/CURRENT_RUN": f"{BASE_RUN_ID}\tPASS\n".encode(),
    }
    for relative, content in unchanged_paths.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)

    acceptance_path = root / f"canonical_publication/provenance/runs/{BASE_RUN_ID}/ACCEPTANCE.json"
    acceptance_path.write_text(
        json.dumps(
            {
                "errors": [],
                "final_manifest": f"canonical_publication/provenance/manifests/{BASE_RUN_ID}.final_artifacts.tsv",
                "run_id": BASE_RUN_ID,
                "status": "PASS",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    unchanged_paths[acceptance_path.relative_to(root).as_posix()] = acceptance_path.read_bytes()

    old_readme = b"old publication documentation\n"
    readme = root / "README.md"
    readme.write_text("current publication documentation\n")
    old_code = b"VALUE = 'old'\n"
    code = root / "canonical_publication/pipeline/src/organelle_pipeline/example.py"
    code.parent.mkdir(parents=True, exist_ok=True)
    code.write_text("VALUE = 'current'\n")

    archive_content = b"original_path\tarchived_path\nold.tsv\tarchive_noncanonical/old.tsv\n"
    archive_copy = root / "canonical_publication/provenance/archive/2026-08-17_pre_remediation/manifest.tsv"
    archive_copy.parent.mkdir(parents=True, exist_ok=True)
    archive_copy.write_bytes(archive_content)

    rows = [
        {
            "canonical_path": relative,
            "size_bytes": str(len(content)),
            "sha256": bytes_digest(content),
            "status": "canonical",
        }
        for relative, content in unchanged_paths.items()
    ]
    rows.extend(
        [
            {
                "canonical_path": "README.md",
                "size_bytes": str(len(old_readme)),
                "sha256": bytes_digest(old_readme),
                "status": "canonical",
            },
            {
                "canonical_path": code.relative_to(root).as_posix(),
                "size_bytes": str(len(old_code)),
                "sha256": bytes_digest(old_code),
                "status": "canonical",
            },
            {
                "canonical_path": "archive_noncanonical/2026-08-17_pre_remediation/manifest.tsv",
                "size_bytes": str(len(archive_content)),
                "sha256": bytes_digest(archive_content),
                "status": "archive_audit_manifest",
            },
        ]
    )
    write_tsv(
        root / f"canonical_publication/provenance/manifests/{BASE_RUN_ID}.final_artifacts.tsv",
        sorted(rows, key=lambda row: row["canonical_path"]),
        ("canonical_path", "size_bytes", "sha256", "status"),
    )


def package_row(path: Path, root: Path, *, status: str = "canonical") -> dict[str, str]:
    return {
        "canonical_path": path.relative_to(root).as_posix(),
        "size_bytes": str(path.stat().st_size),
        "sha256": sha256_file(path),
        "status": status,
    }


def test_package_manifest_validates_regular_file_metadata(tmp_path: Path) -> None:
    artifact = tmp_path / "canonical_publication/results/table.tsv"
    artifact.parent.mkdir(parents=True)
    artifact.write_text("value\n1\n")

    assert validate_package_manifest([package_row(artifact, tmp_path)], tmp_path) == 1


@pytest.mark.parametrize(
    "relative",
    ["/absolute.tsv", "../outside.tsv", "canonical_publication/../outside.tsv", "canonical_publication\\outside.tsv"],
)
def test_package_manifest_rejects_unsafe_paths(tmp_path: Path, relative: str) -> None:
    row = {
        "canonical_path": relative,
        "size_bytes": "1",
        "sha256": "0" * 64,
        "status": "canonical",
    }

    with pytest.raises(PublicationPackageError, match="path"):
        validate_package_manifest([row], tmp_path)


def test_package_manifest_rejects_duplicate_paths(tmp_path: Path) -> None:
    artifact = tmp_path / "canonical_publication/result.tsv"
    artifact.parent.mkdir(parents=True)
    artifact.write_text("value\n")
    row = package_row(artifact, tmp_path)

    with pytest.raises(PublicationPackageError, match="Duplicate"):
        validate_package_manifest([row, row], tmp_path)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("size_bytes", "-1", "size"),
        ("size_bytes", "01", "size"),
        ("sha256", "not-a-digest", "SHA-256"),
        ("status", "unknown", "status"),
    ],
)
def test_package_manifest_rejects_invalid_metadata(tmp_path: Path, field: str, value: str, message: str) -> None:
    artifact = tmp_path / "canonical_publication/result.tsv"
    artifact.parent.mkdir(parents=True)
    artifact.write_text("value\n")
    row = package_row(artifact, tmp_path)
    row[field] = value

    with pytest.raises(PublicationPackageError, match=message):
        validate_package_manifest([row], tmp_path)


def test_package_manifest_rejects_schema_drift_and_symlinks(tmp_path: Path) -> None:
    target = tmp_path / "canonical_publication/target.tsv"
    target.parent.mkdir(parents=True)
    target.write_text("value\n")
    row = package_row(target, tmp_path)
    row["unexpected"] = "field"
    with pytest.raises(PublicationPackageError, match="schema"):
        validate_package_manifest([row], tmp_path)

    link = target.with_name("link.tsv")
    link.symlink_to(target.name)
    link_row = package_row(target, tmp_path)
    link_row["canonical_path"] = link.relative_to(tmp_path).as_posix()
    with pytest.raises(PublicationPackageError, match="regular file"):
        validate_package_manifest([link_row], tmp_path)


def test_create_and_verify_package_without_archive_checkout(tmp_path: Path) -> None:
    miniature_publication(tmp_path)

    created = create_publication_package(tmp_path, BASE_RUN_ID, PACKAGE_ID)
    acceptance = json.loads(created.acceptance.read_text())

    assert not (tmp_path / "archive_noncanonical").exists()
    assert acceptance["status_label"] == "PUBLICATION_PACKAGE_PASS"
    assert acceptance["base_run"]["run_id"] == BASE_RUN_ID
    assert acceptance["archive_snapshot"]["tracked_file_count"] == 1717
    assert verify_publication_package(tmp_path)["package_id"] == PACKAGE_ID


def test_create_rejects_deleted_historical_packaging_file(tmp_path: Path) -> None:
    miniature_publication(tmp_path)
    (tmp_path / "README.md").unlink()

    with pytest.raises(PublicationPackageError, match="historical packaging file"):
        create_publication_package(tmp_path, BASE_RUN_ID, PACKAGE_ID)


def test_package_creation_is_idempotent_but_refuses_changed_content(tmp_path: Path) -> None:
    miniature_publication(tmp_path)
    created = create_publication_package(tmp_path, BASE_RUN_ID, PACKAGE_ID)
    original = {path: path.read_bytes() for path in vars(created).values()}

    repeated = create_publication_package(tmp_path, BASE_RUN_ID, PACKAGE_ID)
    assert {path: path.read_bytes() for path in vars(repeated).values()} == original

    (tmp_path / "README.md").write_text("later packaging change\n")
    with pytest.raises(PublicationPackageError, match="Refusing to replace"):
        create_publication_package(tmp_path, BASE_RUN_ID, PACKAGE_ID)
    assert {path: path.read_bytes() for path in vars(created).values()} == original


def test_verify_rejects_scientific_artifact_tampering(tmp_path: Path) -> None:
    miniature_publication(tmp_path)
    create_publication_package(tmp_path, BASE_RUN_ID, PACKAGE_ID)
    (tmp_path / "canonical_publication/results/run/result.tsv").write_text("result\n2\n")

    with pytest.raises(PublicationPackageError, match="mismatch"):
        verify_publication_package(tmp_path)


def test_cli_verifies_current_package(tmp_path: Path) -> None:
    miniature_publication(tmp_path)
    create_publication_package(tmp_path, BASE_RUN_ID, PACKAGE_ID)
    pipeline_root = Path(__file__).parents[1]
    environment = {**os.environ, "PYTHONPATH": str(pipeline_root / "src")}

    completed = subprocess.run(
        [
            sys.executable,
            str(pipeline_root / "scripts/attest_publication_package.py"),
            "verify",
            "--repository-root",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
        check=False,
        env=environment,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout == f"publication package verified: {PACKAGE_ID}\n"


def test_publication_ci_runs_every_release_gate() -> None:
    repository_root = Path(__file__).resolve().parents[3]
    workflow = (repository_root / ".github/workflows/publication-package.yml").read_text()

    for required in (
        "uses: actions/checkout@08c6903cd8c0fde910a37f88322edcfb5dd907a8",
        "uses: mamba-org/setup-micromamba@d7c9bd84e824b79d2af72a2d4196c7f4300d3476",
        "lfs: true",
        "canonical_publication/environment.yml",
        "pytest -q canonical_publication/pipeline/tests",
        "ruff check canonical_publication/pipeline",
        "ruff format --check canonical_publication/pipeline",
        "attest_publication_package.py verify",
        "git diff --check",
    ):
        assert required in workflow


def test_verify_rejects_acceptance_schema_drift_even_with_updated_index(tmp_path: Path) -> None:
    miniature_publication(tmp_path)
    created = create_publication_package(tmp_path, BASE_RUN_ID, PACKAGE_ID)
    acceptance = json.loads(created.acceptance.read_text())
    acceptance["unexpected"] = "field"
    created.acceptance.write_text(json.dumps(acceptance, indent=2, sort_keys=True) + "\n")
    index_lines = created.checksum_index.read_text().splitlines()
    index_lines[1] = f"{sha256_file(created.acceptance)}  {created.acceptance.relative_to(tmp_path).as_posix()}"
    created.checksum_index.write_text("\n".join(index_lines) + "\n")

    with pytest.raises(PublicationPackageError, match="acceptance schema"):
        verify_publication_package(tmp_path)
