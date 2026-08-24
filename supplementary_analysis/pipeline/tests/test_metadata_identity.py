import threading

from dudleya_supplement.documentation import provider_checksum_gate
from dudleya_supplement.identity import MixedAlleleCall, classify_mixed_allele_samples, index_hopping_status, parse_structured_id
from dudleya_supplement.identity_audit import _candidate_cutoff, _msh_path, _revalidate_provider_row, _sketch_samples
from dudleya_supplement.metadata import apply_metadata_policy, derive_populations


def test_duse_is_retained_for_sample_analyses_but_excluded_from_population_inference() -> None:
    rows = [
        {"sample_id": "DU-173", "popcode": "DUSE", "species": "D. setchellii", "population_name": "source group"},
        {"sample_id": "CY_CAS_LP_328_Du-610", "popcode": "CY_CAS", "species": "D. cymosa", "population_name": "Castle Rock"},
    ]
    corrected, changes = apply_metadata_policy(rows)
    assert corrected[0]["popcode"] == "DUSE"
    assert corrected[0]["sample_analysis_eligible"] == "yes"
    assert corrected[0]["population_inference_eligible"] == "no"
    assert changes == [
        {
            "sample_id": "DU-173",
            "old_popcode": "DUSE",
            "new_popcode_or_EXCLUDED": "EXCLUDED",
            "evidence_source": "repository source-record audit",
            "decision_author": "supplementary pipeline policy",
            "decision_date": "2026-08-24",
            "confidence_or_unresolved": "unresolved",
        }
    ]
    assert [row["popcode"] for row in derive_populations(corrected)] == ["CY_CAS"]


def test_structured_identifier_parser_does_not_call_demultiplex_number_an_index() -> None:
    parsed = parse_structured_id("CY_CAS_LP_328_Du-610_S234_L005_R1_001.fastq.gz")
    assert parsed == {"plate_well": "LP_328", "specimen": "Du-610", "demultiplex_sample": "S234", "lane": "L005"}
    assert index_hopping_status(index_sequences=[], demultiplex_metrics=[]) == "untestable"


def test_mixed_allele_screen_uses_support_fraction_and_robust_outlier_rule() -> None:
    calls = {
        "ordinary": [MixedAlleleCall(20, 18, 2)] * 20,
        "background": [MixedAlleleCall(20, 19, 1)] * 20,
        "suspect": [MixedAlleleCall(20, 10, 10)] * 12 + [MixedAlleleCall(20, 20, 0)] * 8,
    }
    results = {row.sample_id: row for row in classify_mixed_allele_samples(calls)}
    assert results["ordinary"].mixed_site_count == 0
    assert results["suspect"].mixed_site_count == 12
    assert results["suspect"].status == "suspected"


def test_declared_missing_provider_record_is_preserved_without_hashing_repo_root(tmp_path) -> None:
    row = {
        "resolved_source_path": "",
        "expected_md5": "a" * 32,
        "observed_md5": "",
        "status": "DECLARED_MISSING",
    }
    result = _revalidate_provider_row(tmp_path, row)
    assert result["supplementary_observed_md5"] == ""
    assert result["supplementary_status"] == "DECLARED_MISSING_NOT_HASHABLE"


def test_identity_sketches_use_bounded_parallel_workers(tmp_path, monkeypatch) -> None:
    barrier = threading.Barrier(4)
    threads = set()

    def fake_sketch(root, output, sample, inputs):
        threads.add(threading.get_ident())
        barrier.wait(timeout=2)

    monkeypatch.setattr("dudleya_supplement.identity_audit._sketch_sample", fake_sketch)
    rows = [{"sample_id": f"S{index}", "r1_paths": f"S{index}.R1.gz", "r2_paths": f"S{index}.R2.gz"} for index in range(4)]
    paths = _sketch_samples(tmp_path, tmp_path / "work", rows, workers=4)
    assert len(paths) == 4
    assert len(threads) == 4


def test_mash_extension_is_appended_without_collapsing_r1_r2_control_names(tmp_path) -> None:
    assert _msh_path(tmp_path / "sample.R1").name == "sample.R1.msh"
    assert _msh_path(tmp_path / "sample.R2").name == "sample.R2.msh"


def test_self_referential_provider_manifest_is_not_misreported_as_a_checksum_failure(tmp_path) -> None:
    manifest = tmp_path / "md5sum.txt"
    manifest.write_text("self reference cannot have a stable self hash\n")
    row = {
        "resolved_source_path": "md5sum.txt",
        "expected_md5": "a" * 32,
        "observed_md5": "b" * 32,
        "status": "UNVERIFIABLE_SELF_REFERENCE",
    }
    result = _revalidate_provider_row(tmp_path, row)
    assert result["supplementary_observed_md5"]
    assert result["supplementary_status"] == "UNVERIFIABLE_SELF_REFERENCE"


def test_mash_candidate_cutoff_requires_similarity_at_least_as_close_as_best_split_control() -> None:
    assert _candidate_cutoff([0.012, 0.014, 0.036], [0.069]) == 0.012


def test_provider_checksum_gate_accepts_declared_exceptions_but_rejects_resolved_failure() -> None:
    accepted = [
        {"supplementary_status": "PASS"},
        {"supplementary_status": "DECLARED_MISSING_NOT_HASHABLE"},
        {"supplementary_status": "UNVERIFIABLE_SELF_REFERENCE"},
    ]
    assert provider_checksum_gate(accepted) == ("PASS", [])
    assert provider_checksum_gate([*accepted, {"provider_name": "bad.fastq.gz", "supplementary_status": "FAIL"}]) == (
        "FAIL",
        ["bad.fastq.gz"],
    )
