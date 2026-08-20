from pathlib import Path

import pytest
from organelle_pipeline.inventory import (
    ACCEPTABLE_SOURCE_VALIDATION_STATUSES,
    classify_provider_md5,
    inventory_tree,
    source_validation_status,
    validate_inventory,
)


def test_provider_md5_self_reference_is_unverifiable_not_a_data_checksum_failure() -> None:
    assert (
        classify_provider_md5(
            expected="0" * 32,
            observed="1" * 32,
            source_is_provider_manifest=True,
            source_inventory_matches=True,
        )
        == "UNVERIFIABLE_SELF_REFERENCE"
    )


def test_provider_md5_mismatch_still_fails_for_reads() -> None:
    assert (
        classify_provider_md5(
            expected="0" * 32,
            observed="1" * 32,
            source_is_provider_manifest=False,
            source_inventory_matches=True,
        )
        == "FAIL_CHECKSUM"
    )


def test_provider_manifest_self_reference_still_fails_when_source_inventory_drifted() -> None:
    assert (
        classify_provider_md5(
            expected="0" * 32,
            observed="1" * 32,
            source_is_provider_manifest=True,
            source_inventory_matches=False,
        )
        == "FAIL_CHECKSUM"
    )


def test_source_status_reports_both_declared_missing_and_provider_metadata_warning() -> None:
    status = source_validation_status(has_failures=False, has_declared_missing=True, has_self_reference_warning=True)

    assert status == "PASS_WITH_DECLARED_MISSING_AND_PROVIDER_METADATA_WARNING"
    assert status in ACCEPTABLE_SOURCE_VALIDATION_STATUSES


def test_any_true_checksum_failure_makes_source_status_fail() -> None:
    assert source_validation_status(True, True, True) == "FAIL"


def test_archive_inventory_preserves_original_path_and_checksum(tmp_path: Path) -> None:
    snapshot = tmp_path / "archive_noncanonical" / "dated" / "snapshot"
    legacy = snapshot / "old" / "result.tsv"
    legacy.parent.mkdir(parents=True)
    legacy.write_text("value\n1\n")

    records = inventory_tree(
        snapshot,
        repository_root=tmp_path,
        tracked_original_paths={"old/result.tsv"},
        reason="pre-remediation artifact",
    )

    assert len(records) == 1
    record = records[0]
    assert record.original_path == "old/result.tsv"
    assert record.archived_path.endswith("snapshot/old/result.tsv")
    assert record.artifact_type == "file"
    assert record.size_bytes == 8
    assert len(record.sha256) == 64
    assert record.git_status == "tracked"


def test_archive_inventory_records_symlink_without_following_it(tmp_path: Path) -> None:
    snapshot = tmp_path / "archive_noncanonical" / "dated" / "snapshot"
    snapshot.mkdir(parents=True)
    target = snapshot / "target.txt"
    target.write_text("target")
    link = snapshot / "link.txt"
    link.symlink_to("target.txt")

    records = inventory_tree(
        snapshot,
        repository_root=tmp_path,
        tracked_original_paths=set(),
        reason="pre-remediation artifact",
    )

    by_name = {record.original_path: record for record in records}
    assert by_name["link.txt"].artifact_type == "symlink"
    assert by_name["link.txt"].size_bytes == len("target.txt")


def test_archive_inventory_validation_rejects_content_drift(tmp_path: Path) -> None:
    snapshot = tmp_path / "snapshot"
    snapshot.mkdir()
    artifact = snapshot / "result.tsv"
    artifact.write_text("original\n")
    records = inventory_tree(snapshot, tmp_path, set(), "legacy")
    rows = [
        {
            "archived_path": record.archived_path,
            "artifact_type": record.artifact_type,
            "size_bytes": str(record.size_bytes),
            "sha256": record.sha256,
        }
        for record in records
    ]

    assert len(validate_inventory(rows, tmp_path)) == 64
    artifact.write_text("modified\n")
    with pytest.raises(ValueError, match="checksum"):
        validate_inventory(rows, tmp_path)


def test_archive_inventory_validation_rejects_duplicate_paths(tmp_path: Path) -> None:
    artifact = tmp_path / "result.tsv"
    artifact.write_text("original\n")
    records = inventory_tree(tmp_path, tmp_path, set(), "legacy")
    record = next(item for item in records if item.archived_path == "result.tsv")
    row = {
        "archived_path": record.archived_path,
        "artifact_type": record.artifact_type,
        "size_bytes": str(record.size_bytes),
        "sha256": record.sha256,
    }

    with pytest.raises(ValueError, match="Duplicate inventory path"):
        validate_inventory([row, row], tmp_path)
