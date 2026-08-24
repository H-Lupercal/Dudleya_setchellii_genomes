"""Final supplementary acceptance logic."""

from __future__ import annotations


def acceptance_checks(
    *,
    canonical_counts: dict[str, int],
    figure_family_count: int,
    all_artifacts_checksummed: bool,
    canonical_unchanged: bool,
    shared_display_count: int,
    restored_identical_tip_count: int,
    rf_representative_count: int,
) -> dict[str, object]:
    checks: dict[str, dict[str, object]] = {
        "canonical_counts": {
            "status": "PASS" if canonical_counts == {"chloroplast": 276, "mitochondria": 271, "shared": 271} else "FAIL",
            "observed": canonical_counts,
        },
        "figure_family_count": {"status": "PASS" if figure_family_count == 6 else "FAIL", "observed": figure_family_count},
        "artifact_checksums": {"status": "PASS" if all_artifacts_checksummed else "FAIL"},
        "canonical_unchanged": {"status": "PASS" if canonical_unchanged else "FAIL"},
        "shared_display_count": {"status": "PASS" if shared_display_count == 271 else "FAIL", "observed": shared_display_count},
        "restored_identical_tip_count": {
            "status": "PASS" if restored_identical_tip_count == 42 else "FAIL",
            "observed": restored_identical_tip_count,
        },
        "rf_representative_count": {
            "status": "PASS" if rf_representative_count == 229 else "FAIL",
            "observed": rf_representative_count,
        },
    }
    return {"status": "PASS" if all(item["status"] == "PASS" for item in checks.values()) else "FAIL", **checks}
