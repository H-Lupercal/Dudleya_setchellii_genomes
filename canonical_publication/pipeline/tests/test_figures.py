import hashlib
import json
import os
import subprocess
import sys
from io import StringIO
from pathlib import Path

import numpy as np
from Bio import Phylo
from organelle_pipeline.figures import (
    TAXON_COLORS,
    composition_counts,
    distance_aware_layout,
    major_haplotype_labels,
    side_label_layout,
    signed_fst_limit,
    support_is_strong,
    unrooted_tree_layout,
)

EXPECTED_TAXA = (
    "D. abramsii ssp. abramsii",
    "D. abramsii ssp. bettinae",
    "D. abramsii ssp. murina",
    "D. cymosa",
    "D. setchellii",
)


def test_taxon_palette_is_fixed_distinct_and_complete() -> None:
    assert tuple(TAXON_COLORS) == EXPECTED_TAXA
    assert len(set(TAXON_COLORS.values())) == 5
    assert all(color.startswith("#") and len(color) == 7 for color in TAXON_COLORS.values())


def test_signed_fst_limit_is_symmetric_and_preserves_negative_values() -> None:
    values = np.asarray([[np.nan, -0.12, 0.04], [-0.12, np.nan, 0.31], [0.04, 0.31, np.nan]])

    limit = signed_fst_limit(values)

    assert limit == 0.31
    assert (-limit, limit) == (-0.31, 0.31)


def test_signed_fst_limit_rejects_matrix_without_finite_estimates() -> None:
    with np.testing.assert_raises_regex(ValueError, "finite"):
        signed_fst_limit(np.full((2, 2), np.nan))


def test_composition_counts_follow_the_canonical_taxon_order() -> None:
    sample_taxa = {
        "s1": "D. setchellii",
        "s2": "D. abramsii ssp. abramsii",
        "s3": "D. setchellii",
    }

    assert composition_counts(("s1", "s2", "s3"), sample_taxa) == (1, 0, 0, 0, 2)


def test_composition_counts_reject_unknown_taxon() -> None:
    with np.testing.assert_raises_regex(ValueError, "Unknown taxon"):
        composition_counts(("s1",), {"s1": "unresolved"})


def test_haplotype_layout_is_deterministic_and_respects_mutational_distance() -> None:
    edges = (("H1", "H2", 1.0), ("H2", "H3", 5.0))

    first = distance_aware_layout(("H1", "H2", "H3"), edges)
    second = distance_aware_layout(("H1", "H2", "H3"), edges)

    assert first == second
    short_edge = np.linalg.norm(np.asarray(first["H1"]) - np.asarray(first["H2"]))
    long_edge = np.linalg.norm(np.asarray(first["H2"]) - np.asarray(first["H3"]))
    assert long_edge > short_edge * 3


def test_haplotype_labels_are_limited_by_a_declared_sample_count() -> None:
    counts = {"H1": 7, "H2": 1, "H3": 5, "H4": 4}

    assert major_haplotype_labels(counts, minimum_count=5) == ("H1", "H3")


def test_side_label_layout_is_deterministic_balanced_and_nonoverlapping() -> None:
    positions = {"H1": (-0.9, 0.4), "H2": (-0.8, -0.2), "H3": (0.7, 0.1), "H4": (0.9, -0.5)}

    first = side_label_layout(positions, ("H1", "H2", "H3", "H4"))
    second = side_label_layout(positions, ("H4", "H3", "H2", "H1"))

    assert first == second
    assert first["H1"][2] == first["H2"][2] == "right"
    assert first["H3"][2] == first["H4"][2] == "left"
    assert len({first[label][1] for label in ("H1", "H2")}) == 2
    assert len({first[label][1] for label in ("H3", "H4")}) == 2


def test_tree_support_highlight_requires_both_declared_thresholds() -> None:
    assert support_is_strong("80/95")
    assert support_is_strong("100/100")
    assert not support_is_strong("79.9/100")
    assert not support_is_strong("100/94.9")
    assert not support_is_strong("100")
    assert not support_is_strong(None)


def test_unrooted_tree_layout_is_deterministic_and_uses_branch_lengths() -> None:
    tree = Phylo.read(StringIO("((s1:1,s2:1):1,s3:5);"), "newick")

    first = unrooted_tree_layout(tree)
    second = unrooted_tree_layout(tree)

    assert first == second
    tip_positions = {tip.name: np.asarray(first[id(tip)]) for tip in tree.get_terminals()}
    assert np.linalg.norm(tip_positions["s3"] - tip_positions["s1"]) > np.linalg.norm(tip_positions["s1"] - tip_positions["s2"])
    assert np.ptp([position[0] for position in tip_positions.values()]) > 0
    assert np.ptp([position[1] for position in tip_positions.values()]) > 0


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _state(root: Path, path: str, label: str, outputs: tuple[Path, ...], status: str = "complete") -> None:
    state_path = root / path
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state_path.write_text(
        json.dumps(
            {
                "status": status,
                "fingerprint": {"digest": f"{label}-digest"},
                "outputs": {output.relative_to(root).as_posix(): hashlib.sha256(output.read_bytes()).hexdigest() for output in outputs},
            }
        )
    )


def _miniature_figure_repository(root: Path, run_id: str) -> None:
    config = _write(
        root / "canonical_publication/config/publication_config.toml",
        "[qc]\neligibility_depth = 5\nminimum_breadth = 0.80\n",
    )
    cp_fasta = _write(root / "canonical_publication/references/selected/chloroplast.fa", ">chloroplast\n" + "A" * 100 + "\n")
    mt_fasta = _write(root / "canonical_publication/references/selected/mitochondria.fa", ">mitochondria\n" + "A" * 120 + "\n")
    cp_annotation = _write(
        root / "canonical_publication/references/annotations/chloroplast.projected.tsv",
        "feature_id\tfeature_type\tgene\tstart_1based\tend_1based\ncp1\tgene\tmatK\t10\t30\n",
    )
    mt_annotation = _write(
        root / "canonical_publication/references/annotations/mitochondria.projected.tsv",
        "feature_id\tfeature_type\tgene\tstart_1based\tend_1based\nmt1\tgene\tcox1\t25\t60\n",
    )
    cp_ir = _write(
        root / "canonical_publication/references/masks/chloroplast_ir_copies.bed",
        "chloroplast\t0\t10\tIRa\nchloroplast\t90\t100\tIRb\n",
    )
    cp_population = _write(
        root / "canonical_publication/references/masks/chloroplast_population_sites.bed",
        "chloroplast\t10\t90\tpopulation_sites\n",
    )
    mt_repeat = _write(
        root / "canonical_publication/references/masks/mitochondria_repeat_mask.bed",
        "mitochondria\t40\t55\trepeat\n",
    )
    mt_high = _write(
        root / f"canonical_publication/references/masks/{run_id}/mitochondria_high_confidence_sites.bed",
        "mitochondria\t0\t40\thigh_confidence_1\nmitochondria\t55\t120\thigh_confidence_2\n",
    )
    reference_inputs = (cp_fasta, mt_fasta, cp_annotation, mt_annotation, cp_ir, cp_population, mt_repeat)
    _state(root, f"canonical_publication/provenance/runs/{run_id}/references.json", "references", reference_inputs)

    samples = (
        ("s1", "P1", EXPECTED_TAXA[0]),
        ("s2", "P1", EXPECTED_TAXA[0]),
        ("s3", "P2", EXPECTED_TAXA[4]),
        ("s4", "P2", EXPECTED_TAXA[4]),
    )
    metadata_text = "sample_id\tpopcode\tspecies\tpopulation_name\n" + "".join(
        f"{sample}\t{population}\t{taxon}\t{population} name\n" for sample, population, taxon in samples
    )
    metadata_paths = []
    for organelle in ("chloroplast", "mitochondria"):
        metadata_paths.append(_write(root / f"canonical_publication/metadata/qc/{run_id}/{organelle}_samples.tsv", metadata_text))
    breadth = _write(
        root / f"canonical_publication/results/qc/{run_id}/sample_breadth.tsv",
        "sample_id\tcp_unique_sites_breadth_dp5\tmt_unique_sites_breadth_dp5\n"
        "s1\t0.95\t0.93\ns2\t0.90\t0.88\ns3\t0.84\t0.82\ns4\t0.78\t0.76\n",
    )
    preprocessing = _write(
        root / f"canonical_publication/results/qc/{run_id}/read_preprocessing_summary.tsv",
        "sample_id\tread_retention\tpassing_q20_rate\tadapter_trimmed_reads\tinput_reads\tduplication_rate\n"
        "s1\t0.90\t0.98\t10\t100\t0.05\ns2\t0.88\t0.97\t12\t100\t0.06\n"
        "s3\t0.86\t0.96\t8\t100\t0.07\ns4\t0.84\t0.95\t5\t100\t0.08\n",
    )
    _state(
        root,
        f"canonical_publication/provenance/runs/{run_id}/qc.json",
        "qc",
        (breadth, preprocessing, *metadata_paths, mt_high),
    )

    for organelle in ("chloroplast", "mitochondria"):
        pca_coordinates = _write(
            root / f"canonical_publication/results/pca/{run_id}/{organelle}.coordinates.tsv",
            "sample_id\tpopcode\tPC1\tPC2\ns1\tP1\t-1.0\t0.1\ns2\tP1\t-0.8\t-0.1\ns3\tP2\t0.8\t0.2\ns4\tP2\t1.0\t-0.2\n",
        )
        pca_variance = _write(
            root / f"canonical_publication/results/pca/{run_id}/{organelle}.variance.tsv",
            "component\texplained_variance_ratio\nPC1\t0.60\nPC2\t0.25\n",
        )
        _state(
            root,
            f"canonical_publication/provenance/runs/{run_id}/pca/{organelle}.json",
            f"pca-{organelle}",
            (pca_coordinates, pca_variance),
        )
        assignments = _write(
            root / f"canonical_publication/results/haplotypes/{run_id}/{organelle}.sample_haplotypes.tsv",
            "sample_id\tpopcode\thaplotype\ns1\tP1\tH1\ns2\tP1\tH1\ns3\tP2\tH2\ns4\tP2\tH2\n",
        )
        haplotypes = _write(
            root / f"canonical_publication/results/haplotypes/{run_id}/{organelle}.haplotypes.tsv",
            "haplotype\tsample_count\tpopulation_counts\tsequence\nH1\t2\tP1:2\tA\nH2\t2\tP2:2\tG\n",
        )
        edges = _write(
            root / f"canonical_publication/results/haplotypes/{run_id}/{organelle}.network_edges.tsv",
            "haplotype_1\thaplotype_2\tmutational_distance\nH1\tH2\t1\n",
        )
        _state(
            root,
            f"canonical_publication/provenance/runs/{run_id}/haplotypes/{organelle}.json",
            f"haplotypes-{organelle}",
            (assignments, haplotypes, edges),
        )
        population_summary = _write(
            root / f"canonical_publication/results/popgen/{run_id}/{organelle}.population_summary.tsv",
            f"organelle\tpopulation\tsample_count\tnucleotide_diversity\n{organelle}\tP1\t2\t0.01\n{organelle}\tP2\t2\t0.02\n",
        )
        fst = _write(
            root / f"canonical_publication/results/popgen/{run_id}/{organelle}.pairwise_hudson_fst.tsv",
            f"organelle\tpopulation_1\tpopulation_2\thudson_fst\n{organelle}\tP1\tP2\t-0.12\n",
        )
        _state(
            root,
            f"canonical_publication/provenance/runs/{run_id}/popgen/{organelle}.json",
            f"popgen-{organelle}",
            (population_summary, fst),
        )
        tree = _write(
            root / f"canonical_publication/results/trees/{run_id}/{organelle}.primary.treefile",
            "((s1:0.1,s2:0.1)98/100:0.2,(s3:0.1,s4:0.1)96/99:0.2);\n",
        )
        _state(
            root,
            f"canonical_publication/provenance/runs/{run_id}/trees/{organelle}.json",
            f"trees-{organelle}",
            (tree,),
        )
        distance_dir = root / f"canonical_publication/results/supplement/{run_id}/pairwise_distances"
        differences = _write(
            distance_dir / f"{organelle}.sample_pairwise_differences.tsv",
            "sample_id\ts1\ts2\ts3\ts4\ns1\t0\t1\t2\t3\ns2\t1\t0\t3\t2\ns3\t2\t3\t0\t1\ns4\t3\t2\t1\t0\n",
        )
        callable_sites = _write(
            distance_dir / f"{organelle}.sample_pairwise_callable_sites.tsv",
            "sample_id\ts1\ts2\ts3\ts4\ns1\t10\t9\t8\t7\ns2\t9\t10\t7\t8\ns3\t8\t7\t10\t9\ns4\t7\t8\t9\t10\n",
        )
        long_form = _write(
            distance_dir / f"{organelle}.sample_pairwise_distances.tsv",
            f"organelle\tsample_1\tsample_2\tdifferences\tsites_compared\tp_distance\n{organelle}\ts1\ts2\t1\t9\t0.111111111111\n",
        )
        _state(
            root,
            f"canonical_publication/provenance/runs/{run_id}/distances/{organelle}.json",
            f"distances-{organelle}",
            (differences, callable_sites, long_form),
        )
        admixture_dir = root / f"canonical_publication/results/supplement/{run_id}/admixture/{organelle}"
        replicate = _write(
            admixture_dir / "replicate_cv.tsv",
            f"organelle\tk\treplicate\tseed\tcv_error\tq_path\n{organelle}\t1\t1\t1\t1.2\tunused\n{organelle}\t2\t1\t2\t0.9\tunused\n",
        )
        k_summary = _write(
            admixture_dir / "k_summary.tsv",
            "k\tmean_cv_error\tselected\tboundary_optimum\n1\t1.2\tno\tno\n2\t0.9\tyes\tyes\n",
        )
        sample_order = _write(
            admixture_dir / "sample_order.tsv",
            "q_row_1based\tsample_id\n1\ts1\n2\ts2\n3\ts3\n4\ts4\n",
        )
        selected_q = _write(
            admixture_dir / "selected_K2.best_cv.Q.tsv",
            "0.9\t0.1\n0.8\t0.2\n0.2\t0.8\n0.1\t0.9\n",
        )
        selected = _write(
            admixture_dir / "selected_solution.tsv",
            "selected_k\treplicate\tseed\tcv_error\tboundary_optimum\tq_path\n"
            f"2\t1\t2\t0.9\tyes\t{selected_q.relative_to(root).as_posix()}\n",
        )
        _state(
            root,
            f"canonical_publication/provenance/runs/{run_id}/admixture/{organelle}.json",
            f"admixture-{organelle}",
            (replicate, k_summary, sample_order, selected_q, selected),
        )
    _state(
        root,
        f"canonical_publication/provenance/runs/{run_id}/tree_reproducibility.json",
        "treecheck",
        (),
        status="PASS",
    )
    assert config.is_file()


def test_renderer_writes_all_formats_manifest_and_provenance(tmp_path: Path) -> None:
    run_id = "miniature-run"
    _miniature_figure_repository(tmp_path, run_id)
    repository_root = Path(__file__).resolve().parents[3]
    script = repository_root / "canonical_publication/pipeline/scripts/render_figures.py"
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(repository_root / "canonical_publication/pipeline/src")
    environment["MPLCONFIGDIR"] = str(tmp_path / "matplotlib")

    subprocess.run(
        [
            sys.executable,
            str(script),
            "--repository-root",
            str(tmp_path),
            "--config",
            str(tmp_path / "canonical_publication/config/publication_config.toml"),
            "--run-id",
            run_id,
        ],
        cwd=repository_root,
        env=environment,
        check=True,
    )

    figure_dir = tmp_path / f"canonical_publication/reports/figures/{run_id}"
    manifest = figure_dir / "figure_manifest.tsv"
    rows = manifest.read_text().splitlines()
    assert len(rows) == 43
    for extension in ("png", "pdf", "svg"):
        outputs = sorted(figure_dir.glob(f"*.{extension}"))
        assert len(outputs) == 14
        assert all(path.stat().st_size > 100 for path in outputs)
    state = json.loads((tmp_path / f"canonical_publication/provenance/runs/{run_id}/figures.json").read_text())
    assert state["status"] == "complete"
    assert len(state["outputs"]) == 43
    assert all(hashlib.sha256((tmp_path / path).read_bytes()).hexdigest() == digest for path, digest in state["outputs"].items())
