from pathlib import Path

import pytest
from dudleya_supplement.configuration import SupplementConfigurationError, validate_config
from dudleya_supplement.paths import SupplementPathError, assert_output_path
from dudleya_supplement.provenance import (
    StaleSupplementError,
    build_fingerprint,
    code_input_hashes,
    filesystem_snapshot,
    validate_immutable_snapshot,
    validate_resume,
)
from dudleya_supplement.stages import _canonical_fingerprint_value


def valid_config() -> dict[str, object]:
    return {
        "workflow": {"kind": "supplementary", "base_run_id": "publication-20260817", "decision_plan_version": "2.5"},
        "scenarios": {
            "canonical": {"dp": 5, "gq": 20, "missing": 0.20, "breadth": 0.80, "eligibility_dp": 5},
            "permissive": {"dp": 3, "gq": 15, "missing": 0.30, "breadth": 0.70, "eligibility_dp": 3},
            "strict": {"dp": 10, "gq": 30, "missing": 0.10, "breadth": 0.90, "eligibility_dp": 10},
            "mtmask70": {"mask_support": 0.70},
            "mtmask90": {"mask_support": 0.90},
        },
        "seeds": {
            "cp_tree": 271828,
            "mt_tree": 314159,
            "site_resampling": 424200,
            "pi_resampling": 424201,
            "protest": [424210, 424211, 424212, 424213, 424214, 424215],
            "technical_confounders_start": 424300,
        },
        "likelihood_mapping": {"quartets": 100000, "center_limit": 0.15, "side_trigger": 0.20, "split_trigger": 0.20},
        "figures": {
            "families": [
                "robustness",
                "phylogenetic_information",
                "technical_confounders",
                "organelle_comparison",
                "population_diversity",
                "genome_coordinate",
            ]
        },
    }


def test_output_paths_are_confined_to_supplement(tmp_path: Path) -> None:
    root = tmp_path
    supplement = root / "supplementary_analysis/results/table.tsv"
    assert assert_output_path(supplement, root) == supplement.resolve()
    for forbidden in ("canonical_publication/x", "archive_noncanonical/x", "source_data/x", "elsewhere/x"):
        with pytest.raises(SupplementPathError):
            assert_output_path(root / forbidden, root)


def test_configuration_locks_approved_scenarios_and_six_figures() -> None:
    validate_config(valid_config())
    broken = valid_config()
    broken["scenarios"]["strict"]["dp"] = 8  # type: ignore[index]
    with pytest.raises(SupplementConfigurationError, match="strict.dp"):
        validate_config(broken)
    broken = valid_config()
    broken["figures"]["families"] = ["one"]  # type: ignore[index]
    with pytest.raises(SupplementConfigurationError, match="six figure"):
        validate_config(broken)


def test_resume_rejects_changed_inputs() -> None:
    saved = build_fingerprint("metadata", {"input": "aaa"}, {}, ["run metadata"], "commit-a")
    current = build_fingerprint("metadata", {"input": "bbb"}, {}, ["run metadata"], "commit-a")
    with pytest.raises(StaleSupplementError, match="metadata"):
        validate_resume(saved.digest, current)


def test_canonical_snapshot_detects_content_and_timestamp_changes(tmp_path: Path) -> None:
    canonical = tmp_path / "canonical_publication"
    canonical.mkdir()
    file = canonical / "input.tsv"
    file.write_text("a\n")
    saved = filesystem_snapshot(canonical)
    validate_immutable_snapshot(saved, filesystem_snapshot(canonical))
    file.write_text("b\n")
    with pytest.raises(StaleSupplementError, match="immutable canonical"):
        validate_immutable_snapshot(saved, filesystem_snapshot(canonical))


def test_code_fingerprint_includes_supplement_and_explicit_canonical_imports(tmp_path: Path) -> None:
    supplement = tmp_path / "supplementary_analysis/pipeline/src/pkg"
    canonical = tmp_path / "canonical_publication/pipeline/src/base"
    supplement.mkdir(parents=True)
    canonical.mkdir(parents=True)
    (supplement / "a.py").write_text("A = 1\n")
    imported = canonical / "stats.py"
    imported.write_text("B = 2\n")
    hashes = code_input_hashes(tmp_path, [imported])
    assert set(hashes) == {
        "supplementary_analysis/pipeline/src/pkg/a.py",
        "canonical_publication/pipeline/src/base/stats.py",
    }


def test_canonical_fingerprint_reader_accepts_legacy_and_current_state_shapes() -> None:
    digest = "a" * 64
    assert _canonical_fingerprint_value({"fingerprint": digest}) == digest
    assert _canonical_fingerprint_value({"fingerprint": {"digest": digest}}) == digest
    assert _canonical_fingerprint_value({}) == ""
