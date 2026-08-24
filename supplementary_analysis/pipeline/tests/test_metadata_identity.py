from dudleya_supplement.identity import MixedAlleleCall, classify_mixed_allele_samples, index_hopping_status, parse_structured_id
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
