import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from Bio import Phylo
from dudleya_supplement.comparative import expected_pair_count, normalized_unrooted_rf, validate_resampling_spec
from dudleya_supplement.documentation import claim_decision_path
from dudleya_supplement.figures import FIGURE_FAMILIES, validate_figure_manifest
from dudleya_supplement.finalization import resolve_phase2_claims
from dudleya_supplement.rendering import _save
from dudleya_supplement.reporting import acceptance_checks
from organelle_pipeline.popgen import callable_nucleotide_diversity


def test_population_pair_count_formula() -> None:
    assert expected_pair_count(35) == 595
    assert expected_pair_count(1) == 0


def test_resampling_contract_is_exact() -> None:
    validate_resampling_spec(site_draws=1000, site_seed=424200, pi_draws=1000, pi_seed=424201, common_n=4)
    with pytest.raises(ValueError):
        validate_resampling_spec(site_draws=999, site_seed=424200, pi_draws=1000, pi_seed=424201, common_n=4)


def test_figure_manifest_has_exactly_six_families(tmp_path: Path) -> None:
    rows = []
    for index, family in enumerate(FIGURE_FAMILIES, 1):
        for extension in ("png", "pdf", "svg"):
            rows.append({"figure_id": f"S{index}", "family": family, "format": extension, "path": f"x.{extension}"})
    validate_figure_manifest(rows)
    with pytest.raises(ValueError):
        validate_figure_manifest(rows[:-1])


def test_acceptance_requires_canonical_counts_and_checksums() -> None:
    checks = acceptance_checks(
        canonical_counts={"chloroplast": 276, "mitochondria": 271, "shared": 271},
        figure_family_count=6,
        all_artifacts_checksummed=True,
        canonical_unchanged=True,
        shared_display_count=271,
        rf_representative_count=229,
    )
    assert checks["status"] == "PASS"
    checks["canonical_counts"] = {"status": "FAIL"}


def test_callable_denominator_is_recalculated_for_each_sample_draw() -> None:
    first = callable_nucleotide_diversity(["ACGT", "ACNT", "ATGT", "ACGA"])
    second = callable_nucleotide_diversity(["ACGT", "ACNT", "ATGT"])
    assert first.compared_sites != second.compared_sites


def test_multifurcating_unrooted_rf_is_supported(tmp_path: Path) -> None:
    left_path = tmp_path / "left.tree"
    right_path = tmp_path / "right.tree"
    left_path.write_text("(A,B,C,(D,E));\n")
    right_path.write_text("(A,B,(C,D),E);\n")
    left = Phylo.read(left_path, "newick")
    right = Phylo.read(right_path, "newick")
    numerator, denominator, normalized = normalized_unrooted_rf(left, right)
    assert numerator >= 0
    assert denominator >= numerator
    assert 0 <= normalized <= 1


def test_all_figure_formats_use_supported_metadata(tmp_path: Path) -> None:
    figure, axis = plt.subplots()
    axis.plot([0, 1], [0, 1])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        outputs = _save(figure, tmp_path / "figure")
    assert {path.suffix for path in outputs} == {".png", ".pdf", ".svg"}
    assert not [warning for warning in caught if "Unknown infodict keyword" in str(warning.message)]


def test_phase2_claims_are_resolved_from_completed_results() -> None:
    claims = [
        {"metric": "seven_region_likelihood_mapping", "result_status": "PENDING_PHASE2", "required_interpretation_change": ""},
        {"metric": "supported_topology_compatibility", "result_status": "PENDING_PHASE2", "required_interpretation_change": ""},
        {"metric": "resampling_distributions", "result_status": "PENDING_PHASE2", "required_interpretation_change": ""},
    ]
    resolved = resolve_phase2_claims(
        claims,
        likelihood_rows=[
            {
                "organelle": "chloroplast",
                "decision": "TREE_LIKE_NO_NETWORK",
                "center_fraction": "0.013",
                "composition_failed_count": "1",
            },
            {
                "organelle": "mitochondria",
                "decision": "TREE_LIKE_NO_NETWORK",
                "center_fraction": "0.1446",
                "composition_failed_count": "269",
            },
        ],
        rf_row={"rf_numerator": "64", "rf_denominator": "96", "normalized_unrooted_rf": "0.666666666667"},
        resampling_summary={
            "marker_count_result": "observed_within_cp_distribution",
            "sample_size_result": "named_outlier_medians_remain_top_ranked",
        },
    )
    by_metric = {row["metric"]: row for row in resolved}
    assert by_metric["seven_region_likelihood_mapping"]["result_status"] == "PASS_WITH_CAVEAT"
    assert by_metric["supported_topology_compatibility"]["result_status"] == "PASS_WITH_CAVEAT"
    assert by_metric["resampling_distributions"]["result_status"] == "PASS_WITH_CAVEAT"
    assert all(row["result_status"] != "PENDING_PHASE2" for row in resolved)


def test_phase1_and_final_claim_documents_have_separate_owners(tmp_path: Path) -> None:
    phase1 = claim_decision_path(tmp_path, "run", phase1=True)
    final = claim_decision_path(tmp_path, "run", phase1=False)
    assert phase1.name == "claim_analysis_decisions.phase1.tsv"
    assert final.name == "claim_analysis_decisions.tsv"
    assert phase1 != final
