"""Approved comparative, resampling, and confounder analyses."""

from __future__ import annotations

import math
import subprocess
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
from Bio import Phylo, SeqIO
from organelle_pipeline.ordination import prepare_haploid_pca_matrix
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

from .comparative import normalized_unrooted_rf, supported_contracted_tree, validate_resampling_spec
from .io import read_tsv, write_tsv
from .phylogeny import parse_identical_sequence_map
from .scenario import run_scenario_popgen
from .sensitivity import procrustes_permutation_test
from .technical_sensitivity import EXPECTED_FULLY_CALLED_MARKERS, classify_pca_sensitivity, select_fully_called_markers


def _records(path: Path) -> dict[str, str]:
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}


def _bh_adjust(p_values: list[float]) -> list[float]:
    count = len(p_values)
    order = sorted(range(count), key=p_values.__getitem__)
    adjusted = [1.0] * count
    running = 1.0
    for rank_index in range(count - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[index] * count / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def _permutation_spearman(x: np.ndarray, y: np.ndarray, seed: int, permutations: int = 9999) -> tuple[float, float]:
    observed = float(spearmanr(x, y).statistic)
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(permutations):
        value = float(spearmanr(x, y[rng.permutation(len(y))]).statistic)
        exceed += abs(value) >= abs(observed) - 1e-12
    return observed, (exceed + 1) / (permutations + 1)


def finite_pair_spearman(left: list[float], right: list[float]) -> tuple[float, int]:
    """Calculate Spearman agreement after paired finite-value filtering."""
    finite = [(x, y) for x, y in zip(left, right, strict=True) if math.isfinite(x) and math.isfinite(y)]
    if len(finite) < 2:
        return math.nan, len(finite)
    left_finite, right_finite = zip(*finite, strict=True)
    return float(spearmanr(left_finite, right_finite).statistic), len(finite)


def summarize_population_resampling(
    site_rows: list[dict[str, str]],
    pi_rows: list[dict[str, str]],
    *,
    named_outliers: tuple[str, ...] = ("CY_SIE", "CY_CAS"),
) -> dict[str, object]:
    """Summarize the two distinct Figure-5 controls without conflating them."""
    cp_shared = np.asarray([float(row["cp_multi_population_haplotypes"]) for row in site_rows])
    observed_mt = float(site_rows[0]["observed_mt_multi_population_haplotypes"])
    lower, upper = np.quantile(cp_shared, [0.025, 0.975])
    marker_result = "observed_within_cp_distribution" if lower <= observed_mt <= upper else "observed_outside_cp_distribution"

    grouped: dict[str, list[float]] = defaultdict(list)
    sample_sizes = {int(row["sample_size"]) for row in pi_rows}
    if sample_sizes != {4}:
        raise RuntimeError(f"Population pi resampling must use common n=4, found {sorted(sample_sizes)}")
    for row in pi_rows:
        grouped[row["population"]].append(float(row["nucleotide_diversity"]))
    missing = [population for population in named_outliers if population not in grouped]
    if missing:
        raise RuntimeError(f"Named outlier populations missing from pi resampling: {missing}")
    medians = {population: float(np.median(values)) for population, values in grouped.items()}
    ranked = sorted(medians, key=medians.get, reverse=True)
    sample_size_result = (
        "named_outlier_medians_remain_top_ranked"
        if set(ranked[: len(named_outliers)]) == set(named_outliers)
        else "named_outlier_medians_not_top_ranked"
    )
    summary: dict[str, object] = {
        "site_draws": len(site_rows),
        "site_seed": int(site_rows[0]["seed"]) if "seed" in site_rows[0] else 424200,
        "observed_mt_multi_population_haplotypes": observed_mt,
        "cp_draw_q025_multi_population_haplotypes": float(lower),
        "cp_draw_median_multi_population_haplotypes": float(np.median(cp_shared)),
        "cp_draw_q975_multi_population_haplotypes": float(upper),
        "marker_count_result": marker_result,
        "pi_draws_per_population": len(next(iter(grouped.values()))),
        "pi_seed": int(pi_rows[0]["seed"]) if "seed" in pi_rows[0] else 424201,
        "common_sample_size": 4,
        "sample_size_result": sample_size_result,
    }
    for population in named_outliers:
        other_values = [value for other, values in grouped.items() if other != population for value in values]
        summary[f"{population}_median_pi"] = medians[population]
        summary[f"{population}_fraction_other_draws_below_median"] = float(np.mean(np.asarray(other_values) < medians[population]))
    return summary


def run_technical_confounders(root: Path, run_id: str) -> list[Path]:
    output = root / f"supplementary_analysis/results/comparative/{run_id}/technical_confounders.tsv"
    rows: list[dict[str, object]] = []
    seed = 424300
    breadth_rows = {
        row["sample_id"]: row for row in read_tsv(root / "canonical_publication/results/qc/publication-20260817/sample_breadth.tsv")
    }
    for organelle, prefix in (("chloroplast", "cp"), ("mitochondria", "mt")):
        coordinates = read_tsv(root / f"canonical_publication/results/pca/publication-20260817/{organelle}.coordinates.tsv")
        summaries = {
            row["sample_id"]: row
            for row in read_tsv(root / f"canonical_publication/results/alignments/publication-20260817/{organelle}.callable_summary.tsv")
        }
        p_values: list[float] = []
        start = len(rows)
        for pc in ("PC1", "PC2", "PC3"):
            samples = [row["sample_id"] for row in coordinates]
            pc_values = np.asarray([float(row[pc]) for row in coordinates])
            variables = {
                "missingness": np.asarray([1.0 - float(summaries[sample]["callable_fraction_of_analysis_mask"]) for sample in samples]),
                "log_depth": np.asarray(
                    [math.log1p(float(breadth_rows[sample][f"{prefix}_unique_sites_mean_depth"])) for sample in samples]
                ),
                "reference_concordance": np.asarray([float(summaries[sample]["callable_reference_identity"]) for sample in samples]),
            }
            for variable, values in variables.items():
                rho, p_value = _permutation_spearman(pc_values, values, seed)
                rows.append(
                    {
                        "organelle": organelle,
                        "component": pc,
                        "technical_variable": variable,
                        "spearman_rho": f"{rho:.12g}",
                        "permutation_p": f"{p_value:.12g}",
                        "permutations": 9999,
                        "seed": seed,
                        "bh_adjusted_p_within_organelle": "",
                        "nuclear_decoy": "unavailable",
                        "residual_ambiguity": "NUMT/NUPT ambiguity cannot be excluded",
                    }
                )
                p_values.append(p_value)
                seed += 1
        for row, adjusted in zip(rows[start:], _bh_adjust(p_values), strict=True):
            row["bh_adjusted_p_within_organelle"] = f"{adjusted:.12g}"
    write_tsv(output, rows, list(rows[0]), root)
    return [output]


def _vcf_genotypes(path: Path) -> tuple[list[str], np.ndarray]:
    samples = subprocess.run(["bcftools", "query", "-l", str(path)], capture_output=True, text=True, check=True).stdout.splitlines()
    query = subprocess.run(
        ["bcftools", "query", "-f", "%POS[\t%GT]\n", str(path)], capture_output=True, text=True, check=True
    ).stdout.splitlines()
    markers = [[float(value) if value in {"0", "1"} else np.nan for value in line.split("\t")[1:]] for line in query]
    if not markers:
        raise RuntimeError(f"No markers available for fully called-site PCA: {path}")
    return samples, np.asarray(markers, dtype=float).T


def run_complete_site_pca_sensitivity(root: Path, run_id: str, config: dict[str, object]) -> list[Path]:
    output_dir = root / f"supplementary_analysis/results/comparative/{run_id}/technical_sensitivity"
    output_dir.mkdir(parents=True, exist_ok=True)
    protest_seeds = [int(value) for value in config["seeds"]["complete_pca_protest"]]  # type: ignore[index]
    confounder_seed = int(config["seeds"]["complete_pca_confounders_start"])  # type: ignore[index]
    breadth_rows = {
        row["sample_id"]: row for row in read_tsv(root / "canonical_publication/results/qc/publication-20260817/sample_breadth.tsv")
    }
    summary_rows: list[dict[str, object]] = []
    association_rows: list[dict[str, object]] = []
    outputs: list[Path] = []
    for organelle_index, (organelle, prefix) in enumerate((("chloroplast", "cp"), ("mitochondria", "mt"))):
        source = root / f"canonical_publication/results/variants/publication-20260817/{organelle}.mac2_ordination.vcf.gz"
        filtered = output_dir / f"{organelle}.complete_sites.vcf.gz"
        subprocess.run(["bcftools", "view", "-i", "F_MISSING=0", "-Oz", "-o", str(filtered), str(source)], check=True)
        subprocess.run(["bcftools", "index", "-f", str(filtered)], check=True)
        samples, genotypes = _vcf_genotypes(filtered)
        complete = select_fully_called_markers(genotypes)
        if complete.shape[1] != EXPECTED_FULLY_CALLED_MARKERS[organelle]:
            raise RuntimeError(
                f"Expected {EXPECTED_FULLY_CALLED_MARKERS[organelle]} fully called {organelle} MAC>=2 SNPs, found {complete.shape[1]}"
            )
        matrix = prepare_haploid_pca_matrix(complete)
        component_count = min(10, matrix.shape[0] - 1, matrix.shape[1])
        model = PCA(n_components=component_count, svd_solver="full")
        coordinates = model.fit_transform(matrix)
        canonical_rows = {
            row["sample_id"]: row
            for row in read_tsv(root / f"canonical_publication/results/pca/publication-20260817/{organelle}.coordinates.tsv")
        }
        if set(samples) != set(canonical_rows):
            raise RuntimeError(f"Fully called-site PCA sample mismatch for {organelle}")
        coordinate_path = output_dir / f"{organelle}.complete_sites.coordinates.tsv"
        variance_path = output_dir / f"{organelle}.complete_sites.variance.tsv"
        coordinate_rows = [
            {
                "sample_id": sample,
                "popcode": canonical_rows[sample]["popcode"],
                **{f"PC{index + 1}": f"{value:.12g}" for index, value in enumerate(values)},
            }
            for sample, values in zip(samples, coordinates, strict=True)
        ]
        write_tsv(
            coordinate_path,
            coordinate_rows,
            ["sample_id", "popcode", *[f"PC{index}" for index in range(1, component_count + 1)]],
            root,
        )
        write_tsv(
            variance_path,
            [
                {"component": f"PC{index}", "explained_variance_ratio": f"{value:.12g}"}
                for index, value in enumerate(model.explained_variance_ratio_, 1)
            ],
            ["component", "explained_variance_ratio"],
            root,
        )
        canonical_matrix = np.asarray(
            [[float(canonical_rows[sample][f"PC{component}"]) for component in range(1, 4)] for sample in samples]
        )
        protest = procrustes_permutation_test(
            canonical_matrix,
            coordinates[:, :3],
            permutations=9999,
            seed=protest_seeds[organelle_index],
        )
        status = classify_pca_sensitivity(protest.correlation, protest.p_value)
        summary_rows.append(
            {
                "organelle": organelle,
                "samples": len(samples),
                "canonical_mac2_markers": int(
                    subprocess.run(["bcftools", "view", "-H", str(source)], capture_output=True, text=True, check=True).stdout.count("\n")
                ),
                "fully_called_mac2_markers": complete.shape[1],
                "protest_r": f"{protest.correlation:.12g}",
                "protest_p": f"{protest.p_value:.12g}",
                "permutations": protest.permutations,
                "seed": protest.seed,
                "status": status,
            }
        )
        summaries = {
            row["sample_id"]: row
            for row in read_tsv(root / f"canonical_publication/results/alignments/publication-20260817/{organelle}.callable_summary.tsv")
        }
        start = len(association_rows)
        p_values: list[float] = []
        for component in range(1, 4):
            pc_values = coordinates[:, component - 1]
            variables = {
                "missingness": np.asarray([1.0 - float(summaries[sample]["callable_fraction_of_analysis_mask"]) for sample in samples]),
                "log_depth": np.asarray(
                    [math.log1p(float(breadth_rows[sample][f"{prefix}_unique_sites_mean_depth"])) for sample in samples]
                ),
                "reference_concordance": np.asarray([float(summaries[sample]["callable_reference_identity"]) for sample in samples]),
            }
            for variable, values in variables.items():
                rho, p_value = _permutation_spearman(pc_values, values, confounder_seed)
                association_rows.append(
                    {
                        "organelle": organelle,
                        "component": f"PC{component}",
                        "technical_variable": variable,
                        "spearman_rho": f"{rho:.12g}",
                        "permutation_p": f"{p_value:.12g}",
                        "permutations": 9999,
                        "seed": confounder_seed,
                        "bh_adjusted_p_within_organelle": "",
                        "interpretation_limit": ("association cannot distinguish genuine divergence from reference-mapping bias"),
                    }
                )
                p_values.append(p_value)
                confounder_seed += 1
        for row, adjusted in zip(association_rows[start:], _bh_adjust(p_values), strict=True):
            row["bh_adjusted_p_within_organelle"] = f"{adjusted:.12g}"
        outputs.extend([filtered, Path(f"{filtered}.csi"), coordinate_path, variance_path])
    summary_path = output_dir / "complete_site_pca_summary.tsv"
    associations_path = output_dir / "complete_site_technical_confounders.tsv"
    write_tsv(summary_path, summary_rows, list(summary_rows[0]), root)
    write_tsv(associations_path, association_rows, list(association_rows[0]), root)
    return [*outputs, summary_path, associations_path]


def _representative_names(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    expected = int(lines[0].split()[0])
    names = [line.split()[0] for line in lines[1:] if line.strip()]
    if len(names) != expected:
        raise RuntimeError(f"Expected {expected} representative names, found {len(names)}")
    return names


def run_organelle_comparison(root: Path, run_id: str) -> list[Path]:
    output_dir = root / f"supplementary_analysis/results/comparative/{run_id}/organelle_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)
    representatives = set(
        _representative_names(root / "canonical_publication/results/trees/publication-20260817/mitochondria.primary.uniqueseq.phy")
    )
    if len(representatives) != 229:
        raise RuntimeError(f"Expected 229 mitochondrial representatives, found {len(representatives)}")
    cp_tree = supported_contracted_tree(
        root / "canonical_publication/results/trees/publication-20260817/chloroplast.primary.treefile", representatives
    )
    mt_tree = supported_contracted_tree(
        root / "canonical_publication/results/trees/publication-20260817/mitochondria.primary.treefile", representatives
    )
    rf_numerator, rf_denominator, normalized = normalized_unrooted_rf(cp_tree, mt_tree)
    cp_out = output_dir / "chloroplast.support_contracted_229.treefile"
    mt_out = output_dir / "mitochondria.support_contracted_229.treefile"
    Phylo.write(cp_tree, cp_out, "newick")
    Phylo.write(mt_tree, mt_out, "newick")
    rf_path = output_dir / "supported_unrooted_rf.tsv"
    write_tsv(
        rf_path,
        [
            {
                "taxon_space": "229_mitochondrial_unique_representatives",
                "rf_numerator": rf_numerator,
                "rf_denominator": rf_denominator,
                "normalized_unrooted_rf": f"{normalized:.12g}",
                "branch_rule": "retain only SH-aLRT>=80 and UFBoot>=95",
                "interpretation": "supported-topology compatibility; not total evolutionary disagreement",
            }
        ],
        ["taxon_space", "rf_numerator", "rf_denominator", "normalized_unrooted_rf", "branch_rule", "interpretation"],
        root,
    )
    cp_full = Phylo.read(root / "canonical_publication/results/trees/publication-20260817/chloroplast.primary.treefile", "newick")
    mt_full = Phylo.read(root / "canonical_publication/results/trees/publication-20260817/mitochondria.primary.treefile", "newick")
    cp_order = {tip.name: index for index, tip in enumerate(cp_full.get_terminals(), 1)}
    mt_order = {tip.name: index for index, tip in enumerate(mt_full.get_terminals(), 1)}
    shared = sorted(set(cp_order) & set(mt_order))
    mapping = parse_identical_sequence_map(
        (root / "canonical_publication/results/trees/publication-20260817/mitochondria.primary.log").read_text()
    )
    tangle = output_dir / "tanglegram_271_tip_mapping.tsv"
    write_tsv(
        tangle,
        [
            {
                "sample_id": sample,
                "chloroplast_tip_order": cp_order[sample],
                "mitochondria_tip_order": mt_order[sample],
                "mitochondria_representative": mapping.get(sample, sample),
                "identical_zero_length_tip_group": "yes" if sample in mapping else "no",
            }
            for sample in shared
        ],
        ["sample_id", "chloroplast_tip_order", "mitochondria_tip_order", "mitochondria_representative", "identical_zero_length_tip_group"],
        root,
    )
    if len(shared) != 271:
        raise RuntimeError(f"Tanglegram must display 271 shared samples, found {len(shared)}")
    fst: dict[str, dict[tuple[str, str], float]] = {}
    for organelle in ("chloroplast", "mitochondria"):
        values = {}
        for row in read_tsv(
            root / f"supplementary_analysis/results/sensitivity/{run_id}/canonical/popgen/{organelle}.pairwise_hudson_fst.tsv"
        ):
            values[(row["population_1"], row["population_2"])] = float(row["hudson_fst"])
        fst[organelle] = values
    pairs = sorted(set(fst["chloroplast"]) & set(fst["mitochondria"]))
    cp_values = [fst["chloroplast"][pair] for pair in pairs]
    mt_values = [fst["mitochondria"][pair] for pair in pairs]
    rho, finite_count = finite_pair_spearman(cp_values, mt_values)
    fst_path = output_dir / "common_pair_fst_agreement.tsv"
    write_tsv(
        fst_path,
        [
            {
                "population_1": pair[0],
                "population_2": pair[1],
                "chloroplast_hudson_fst": f"{fst['chloroplast'][pair]:.12g}",
                "mitochondria_hudson_fst": f"{fst['mitochondria'][pair]:.12g}",
                "global_spearman_rho": f"{rho:.12g}",
                "eligible_common_pair_count": len(pairs),
                "finite_common_pair_count": finite_count,
                "nonfinite_common_pair_count": len(pairs) - finite_count,
                "used_for_correlation": (
                    "yes" if math.isfinite(fst["chloroplast"][pair]) and math.isfinite(fst["mitochondria"][pair]) else "no"
                ),
            }
            for pair in pairs
        ],
        [
            "population_1",
            "population_2",
            "chloroplast_hudson_fst",
            "mitochondria_hudson_fst",
            "global_spearman_rho",
            "eligible_common_pair_count",
            "finite_common_pair_count",
            "nonfinite_common_pair_count",
            "used_for_correlation",
        ],
        root,
    )
    return [cp_out, mt_out, rf_path, tangle, fst_path]


def _haplotype_sharing(sequences: dict[str, str], metadata: dict[str, str]) -> tuple[int, int, int]:
    assigned = {sample: sequence for sample, sequence in sequences.items() if "N" not in sequence and "-" not in sequence}
    populations: dict[str, set[str]] = defaultdict(set)
    for sample, sequence in assigned.items():
        populations[sequence].add(metadata[sample])
    return len(assigned), len(populations), sum(len(values) > 1 for values in populations.values())


def run_population_resampling(root: Path, run_id: str) -> list[Path]:
    validate_resampling_spec(site_draws=1000, site_seed=424200, pi_draws=1000, pi_seed=424201, common_n=4)
    output_dir = root / f"supplementary_analysis/results/comparative/{run_id}/population_diversity"
    output_dir.mkdir(parents=True, exist_ok=True)
    cp = _records(root / "canonical_publication/results/alignments/publication-20260817/chloroplast.callable_alignment.fa")
    mt = _records(root / "canonical_publication/results/alignments/publication-20260817/mitochondria.callable_alignment.fa")
    metadata_rows = read_tsv(root / "supplementary_analysis/metadata/samples/samples.corrected-20260824.tsv")
    metadata = {row["sample_id"]: row["popcode"] for row in metadata_rows}
    shared = sorted(set(cp) & set(mt))
    cp_positions = [
        int(value) - 1
        for value in subprocess.run(
            [
                "bcftools",
                "query",
                "-f",
                "%POS\\n",
                str(root / "canonical_publication/results/variants/publication-20260817/chloroplast.primary.vcf.gz"),
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.splitlines()
    ]
    if len(cp_positions) < 146:
        raise RuntimeError("Chloroplast has fewer than the required 146 primary sites")
    mt_positions = [
        int(value) - 1
        for value in subprocess.run(
            [
                "bcftools",
                "query",
                "-f",
                "%POS\\n",
                str(root / "canonical_publication/results/variants/publication-20260817/mitochondria.primary.vcf.gz"),
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.splitlines()
    ]
    mt_subset = {sample: "".join(mt[sample][position] for position in mt_positions) for sample in shared}
    mt_assigned, mt_haplotypes, mt_shared = _haplotype_sharing(mt_subset, metadata)
    rng = np.random.default_rng(424200)
    site_rows = []
    for draw in range(1, 1001):
        positions = rng.choice(cp_positions, size=146, replace=False)
        subset = {sample: "".join(cp[sample][position] for position in positions) for sample in shared}
        assigned, haplotypes, multi = _haplotype_sharing(subset, metadata)
        site_rows.append(
            {
                "draw": draw,
                "seed": 424200,
                "chloroplast_sites": 146,
                "cp_assigned_samples": assigned,
                "cp_haplotype_count": haplotypes,
                "cp_multi_population_haplotypes": multi,
                "observed_mt_assigned_samples": mt_assigned,
                "observed_mt_haplotype_count": mt_haplotypes,
                "observed_mt_multi_population_haplotypes": mt_shared,
                "interpretation": "marker-count sensitivity only",
            }
        )
    site_path = output_dir / "chloroplast_146_site_resampling.tsv"
    write_tsv(site_path, site_rows, list(site_rows[0]), root)
    corrected = {row["sample_id"]: row for row in metadata_rows}
    population_samples: dict[str, list[str]] = defaultdict(list)
    for sample in cp:
        if corrected[sample]["population_inference_eligible"] == "yes":
            population_samples[metadata[sample]].append(sample)
    pair_stats: dict[tuple[str, str], tuple[int, int]] = {}
    for samples in population_samples.values():
        for left, right in combinations(samples, 2):
            callable_sites = differences = 0
            for a, b in zip(cp[left], cp[right], strict=True):
                if a in "ACGT" and b in "ACGT":
                    callable_sites += 1
                    differences += a != b
            pair_stats[tuple(sorted((left, right)))] = (differences, callable_sites)
    rng = np.random.default_rng(424201)
    pi_rows = []
    for population, samples in sorted(population_samples.items()):
        if len(samples) < 4:
            continue
        for draw in range(1, 1001):
            chosen = sorted(rng.choice(samples, size=4, replace=False).tolist())
            values = [pair_stats[tuple(sorted(pair))] for pair in combinations(chosen, 2)]
            differences = sum(value[0] for value in values)
            denominator = sum(value[1] for value in values)
            pi_rows.append(
                {
                    "population": population,
                    "draw": draw,
                    "seed": 424201,
                    "sample_size": 4,
                    "sample_ids": ",".join(chosen),
                    "pairwise_differences": differences,
                    "pairwise_callable_sites": denominator,
                    "nucleotide_diversity": f"{differences / denominator:.12g}" if denominator else "nan",
                }
            )
    pi_path = output_dir / "population_pi_n4_resampling.tsv"
    write_tsv(pi_path, pi_rows, list(pi_rows[0]), root)
    summary = summarize_population_resampling(site_rows, pi_rows)
    summary_path = output_dir / "population_resampling_summary.tsv"
    write_tsv(summary_path, [summary], list(summary), root)
    return [site_path, pi_path, summary_path]


def _depth_windows(root: Path, run_id: str, organelle: str, samples: list[str], length: int, window: int) -> np.ndarray:
    cache = root / f"supplementary_analysis/work/{run_id}/coordinate_depth/{organelle}.mean_depth.npy"
    if cache.is_file():
        return np.load(cache)
    cache.parent.mkdir(parents=True, exist_ok=True)
    totals = np.zeros(length, dtype=np.float64)
    for sample in samples:
        bam = root / f"canonical_publication/work/publication-20260817/mapping/{sample}.organelle.bam"
        process = subprocess.Popen(
            ["samtools", "depth", "-aa", "-q", "20", "-Q", "20", "-r", organelle, str(bam)],
            stdout=subprocess.PIPE,
            text=True,
        )
        assert process.stdout is not None
        for line in process.stdout:
            _, position, depth = line.rstrip().split("\t")
            totals[int(position) - 1] += int(depth)
        if process.wait() != 0:
            raise RuntimeError(f"samtools depth failed for {sample} {organelle}")
    totals /= len(samples)
    np.save(cache, totals)
    return totals


def _bed_mask(path: Path, record: str, length: int) -> np.ndarray:
    values = np.zeros(length, dtype=bool)
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        name, start, end, *_ = line.split("\t")
        if name == record:
            values[int(start) : int(end)] = True
    return values


def run_coordinate_tracks(root: Path, run_id: str) -> list[Path]:
    output = root / f"supplementary_analysis/results/comparative/{run_id}/genome_coordinate_windows.tsv"
    rows: list[dict[str, object]] = []
    for organelle, length in (("chloroplast", 150274), ("mitochondria", 243359)):
        sample_rows = read_tsv(root / f"canonical_publication/metadata/qc/publication-20260817/{organelle}_samples.tsv")
        samples = [row["sample_id"] for row in sample_rows]
        mask_path = (
            root / "canonical_publication/references/masks/chloroplast_population_sites.bed"
            if organelle == "chloroplast"
            else root / "canonical_publication/references/masks/publication-20260817/mitochondria_high_confidence_sites.bed"
        )
        callable_mask = _bed_mask(mask_path, organelle, length)
        repeat_path = (
            root / "canonical_publication/references/masks/chloroplast_duplicate_ir_mask.bed"
            if organelle == "chloroplast"
            else root / "canonical_publication/references/masks/mitochondria_repeat_mask.bed"
        )
        repeat_mask = _bed_mask(repeat_path, organelle, length)
        variants = set(
            int(value) - 1
            for value in subprocess.run(
                [
                    "bcftools",
                    "query",
                    "-f",
                    "%POS\\n",
                    str(root / f"canonical_publication/results/variants/publication-20260817/{organelle}.primary.vcf.gz"),
                ],
                capture_output=True,
                text=True,
                check=True,
            ).stdout.splitlines()
        )
        mean_depth = _depth_windows(root, run_id, organelle, samples, length, 5000)
        for start in range(0, length, 5000):
            end = min(length, start + 5000)
            callable_bases = int(callable_mask[start:end].sum())
            count = sum(start <= position < end for position in variants)
            rows.append(
                {
                    "organelle": organelle,
                    "start_0based": start,
                    "end_0based_exclusive": end,
                    "window_bp": end - start,
                    "callable_bases": callable_bases,
                    "repeat_or_duplicate_ir_bases": int(repeat_mask[start:end].sum()),
                    "primary_snp_count": count,
                    "variation_per_callable_kb": f"{count / (callable_bases / 1000):.12g}" if callable_bases else "nan",
                    "mean_filtered_depth": f"{float(mean_depth[start:end].mean()):.12g}",
                    "non_overlapping_window": "yes",
                }
            )
    write_tsv(output, rows, list(rows[0]), root)
    return [output]


def run_comparative_analyses(root: Path, run_id: str, config: dict[str, object]) -> list[Path]:
    return [
        *run_identity_sensitivity(root, run_id),
        *run_technical_confounders(root, run_id),
        *run_complete_site_pca_sensitivity(root, run_id, config),
        *run_organelle_comparison(root, run_id),
        *run_population_resampling(root, run_id),
        *run_coordinate_tracks(root, run_id),
    ]


def run_identity_sensitivity(root: Path, run_id: str) -> list[Path]:
    """Recalculate population summaries after excluding screen-suspected samples."""
    outcomes = {
        row["sample_id"]: row["outcome"]
        for row in read_tsv(root / f"supplementary_analysis/results/verification/{run_id}/identity/sample_identity_outcomes.tsv")
    }
    metadata = read_tsv(root / "supplementary_analysis/metadata/samples/samples.corrected-20260824.tsv")
    for row in metadata:
        if outcomes.get(row["sample_id"]) == "suspected":
            row["population_inference_eligible"] = "no"
            row["population_exclusion_reason"] = "identity screen suspected; with-without sensitivity"
    output_dir = root / f"supplementary_analysis/results/comparative/{run_id}/identity_sensitivity"
    metadata_path = output_dir / "samples.without_suspected.tsv"
    write_tsv(metadata_path, metadata, list(metadata[0]), root)
    outputs = [metadata_path]
    for organelle in ("chloroplast", "mitochondria"):
        eligible = {
            row["sample_id"] for row in read_tsv(root / f"canonical_publication/metadata/qc/publication-20260817/{organelle}_samples.tsv")
        }
        organelle_metadata = [row for row in metadata if row["sample_id"] in eligible]
        organelle_metadata_path = output_dir / f"{organelle}_samples.without_suspected.tsv"
        write_tsv(organelle_metadata_path, organelle_metadata, list(organelle_metadata[0]), root)
        outputs.append(organelle_metadata_path)
        population, fst = run_scenario_popgen(
            root,
            run_id,
            "identity_without_suspected",
            organelle,
            root / f"canonical_publication/results/alignments/publication-20260817/{organelle}.callable_alignment.fa",
            organelle_metadata_path,
            root / f"canonical_publication/results/haplotypes/publication-20260817/{organelle}.sample_haplotypes.tsv",
        )
        outputs.extend([population, fst])
    return outputs
