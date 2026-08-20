from organelle_pipeline.analysis import (
    alignment_callability_counts,
    alignment_site_counts,
    build_iqtree_command,
    is_strong_iqtree_support,
    parse_iqtree_support,
    select_best_k,
)


def test_pca_stage_emits_data_not_figures() -> None:
    root = __import__("pathlib").Path(__file__).resolve().parents[3]
    source = (root / "canonical_publication/pipeline/scripts/run_pca.py").read_text()

    assert "matplotlib" not in source
    assert "figure_path" not in source
    assert "reports/figures" not in source


def test_haplotype_stage_emits_data_not_figures() -> None:
    root = __import__("pathlib").Path(__file__).resolve().parents[3]
    source = (root / "canonical_publication/pipeline/scripts/run_haplotypes.py").read_text()

    assert "matplotlib" not in source
    assert "figure_path" not in source
    assert "reports/figures" not in source


def test_variable_sites_are_not_mislabeled_as_parsimony_informative() -> None:
    records = {
        "s1": "AAA",
        "s2": "ATA",
        "s3": "ATT",
        "s4": "ATT",
    }
    # Position 2 is a singleton variable site; position 3 is 2-vs-2 and
    # therefore parsimony-informative.
    assert alignment_site_counts(records) == (2, 1)


def test_alignment_contribution_separates_coordinate_span_from_callable_sites() -> None:
    summary = alignment_callability_counts(
        {
            "s1": "AN",
            "s2": "AN",
            "s3": "AT",
        }
    )
    assert summary.coordinate_span_sites == 2
    assert summary.sites_with_any_callable_sample == 2
    assert summary.sites_with_at_least_two_callable_samples == 1
    assert summary.jointly_callable_sites == 1


def test_admixture_endpoint_optimum_is_explicit() -> None:
    choice = select_best_k({1: 1.2, 2: 1.1, 12: 0.9}, tested_min=1, tested_max=12)
    assert choice.k == 12
    assert choice.is_boundary


def test_primary_tree_uses_model_finder_standard_support_and_seed() -> None:
    command = build_iqtree_command("cp.fa", "cp", seed=271828)
    assert "-m MFP" in command
    assert "-st DNA" in command
    assert "-alrt 1000" in command
    assert "-B 1000" in command
    assert "-bnni" in command
    assert "-seed 271828" in command
    assert "10000" not in command


def test_iqtree_dual_support_labels_are_not_misread_as_one_bootstrap_value() -> None:
    support = parse_iqtree_support("82.5/97", None)
    assert support.sh_alrt == 82.5
    assert support.ultrafast_bootstrap == 97.0
    assert is_strong_iqtree_support(support)
    assert not is_strong_iqtree_support(support, minimum_sh_alrt=90.0, minimum_ultrafast_bootstrap=95.0)
    assert not is_strong_iqtree_support(parse_iqtree_support("79.9/100", None))
