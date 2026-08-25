"""Deterministic rendering of the six approved figure families."""

from __future__ import annotations

import subprocess
from datetime import UTC, datetime
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .figures import FIGURE_FAMILIES, FORMATS, validate_figure_manifest
from .io import read_tsv, write_tsv

COLORS = {"chloroplast": "#27864a", "mitochondria": "#8a4f9e"}


def _save(fig, prefix: Path) -> list[Path]:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = []
    for extension in FORMATS:
        path = prefix.with_suffix(f".{extension}")
        kwargs = {"dpi": 300} if extension == "png" else {}
        metadata: dict[str, object] = {"Creator": "Dudleya supplementary pipeline"}
        if extension == "pdf":
            metadata["CreationDate"] = datetime(2026, 8, 24, tzinfo=UTC)
        elif extension == "svg":
            metadata["Date"] = "2026-08-24"
        fig.savefig(path, bbox_inches="tight", metadata=metadata, **kwargs)
        outputs.append(path)
    plt.close(fig)
    return outputs


def _robustness(root: Path, run_id: str, prefix: Path) -> list[Path]:
    rows = read_tsv(root / f"supplementary_analysis/results/sensitivity/{run_id}/sensitivity_summary.tsv")
    labels = [f"{row['scenario']}\n{'cp' if row['organelle'] == 'chloroplast' else 'mt'}" for row in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    axes[0, 0].bar(x, [float(row["pi_spearman_rho"]) for row in rows], color="#4c78a8")
    axes[0, 0].axhline(0.95, color="black", linestyle="--", linewidth=1)
    axes[0, 0].set(title="π rank agreement", ylabel="Spearman ρ", ylim=(0, 1.03))
    axes[0, 1].bar(x, [float(row["fst_spearman_rho"]) for row in rows], color="#f58518")
    axes[0, 1].axhline(0.95, color="black", linestyle="--", linewidth=1)
    axes[0, 1].set(title="FST rank agreement", ylabel="Spearman ρ", ylim=(0, 1.03))
    axes[1, 0].bar(x, [float(row["protest_r"]) for row in rows], color="#54a24b")
    axes[1, 0].axhline(0.90, color="black", linestyle="--", linewidth=1)
    axes[1, 0].set(title="PC1–PC3 Procrustes agreement", ylabel="r", ylim=(0, 1.03))
    for row in rows:
        short = "cp" if row["organelle"] == "chloroplast" else "mt"
        axes[1, 1].scatter(
            int(row["eligible_samples"]),
            int(row["retained_snps"]),
            color=COLORS[row["organelle"]],
            s=55,
            label=f"{row['scenario']} {short}",
        )
    axes[1, 1].set(title="Samples and retained SNPs", xlabel="Eligible samples", ylabel="SNPs")
    axes[1, 1].legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    for axis in axes.flat[:3]:
        axis.set_xticks(x, labels, rotation=35, ha="right", fontsize=7)
    fig.suptitle("Supplementary Figure 1 — Filtering robustness")
    fig.tight_layout()
    return _save(fig, prefix)


def _phylogenetic_information(root: Path, run_id: str, prefix: Path) -> list[Path]:
    rows = read_tsv(root / f"supplementary_analysis/results/phylogeny/{run_id}/likelihood_mapping/likelihood_mapping_summary.tsv")
    sensitivity = read_tsv(
        root / f"supplementary_analysis/results/phylogeny/{run_id}/likelihood_mapping_sensitivity/likelihood_mapping_sensitivity.tsv"
    )[0]
    panels = [*rows, sensitivity]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
    palette = ["#2e8b57", "#4aa564", "#72bd72", "#d6a140", "#e4ba66", "#efcf91", "#a7a7a7"]
    for axis, row in zip(axes, panels, strict=True):
        bottom = 0.0
        for index in range(1, 8):
            value = float(row[f"region_{index}_fraction"])
            axis.bar([0], [value], bottom=bottom, color=palette[index - 1], label=f"Region {index}")
            bottom += value
        axis.set(
            xticks=[],
            ylim=(0, 1),
            ylabel="Quartet fraction",
            xlabel=(
                f"Composition failures: {row['composition_failed_count']}/{row['alignment_sequence_count']}; "
                f">50% gaps/ambiguity: {row['over_50pct_ambiguity_count']}"
            ),
        )
        title = row["organelle"].capitalize()
        decision = row.get("decision", "DIAGNOSTIC_ONLY").replace("_", " ").title()
        if row.get("analysis") == "mitochondria_mask_restricted":
            title = "Mitochondria\n43,182-base mask sensitivity"
            decision = "Diagnostic only — no network trigger"
        axis.set_title(title, y=1.10)
        axis.text(0.5, 1.02, decision, transform=axis.transAxes, ha="center", fontsize=9)
    axes[2].legend(ncol=2, fontsize=8, loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.suptitle("Supplementary Figure 2 — Seven-region likelihood mapping\nPrimary results plus diagnostic-only mt mask sensitivity")
    fig.tight_layout()
    return _save(fig, prefix)


def _genotype_matrix(vcf: Path) -> np.ndarray:
    text = subprocess.run(["bcftools", "query", "-f", "%POS[\\t%GT]\\n", str(vcf)], capture_output=True, text=True, check=True).stdout
    rows = []
    for line in text.splitlines():
        rows.append([float(value) if value in {"0", "1"} else np.nan for value in line.split("\t")[1:]])
    return np.asarray(rows, dtype=float).T


def _technical(root: Path, run_id: str, prefix: Path) -> list[Path]:
    stats = read_tsv(root / f"supplementary_analysis/results/comparative/{run_id}/technical_confounders.tsv")
    complete_stats = read_tsv(
        root / f"supplementary_analysis/results/comparative/{run_id}/technical_sensitivity/complete_site_technical_confounders.tsv"
    )
    complete_summary = {
        row["organelle"]: row
        for row in read_tsv(
            root / f"supplementary_analysis/results/comparative/{run_id}/technical_sensitivity/complete_site_pca_summary.tsv"
        )
    }
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for axis, organelle in zip(axes[0], ("chloroplast", "mitochondria"), strict=True):
        matrix = _genotype_matrix(root / f"canonical_publication/results/variants/publication-20260817/{organelle}.primary.vcf.gz")
        show = matrix[:, : min(500, matrix.shape[1])]
        axis.imshow(show, aspect="auto", interpolation="nearest", cmap="viridis", vmin=0, vmax=1)
        axis.set(
            title=f"{organelle.capitalize()} genotypes (first 500 SNPs)",
            xlabel="SNP (purple=reference; yellow=alternate; white=missing)",
            ylabel="Sample",
        )
    for axis, organelle in zip(axes[1], ("chloroplast", "mitochondria"), strict=True):
        selected = [row for row in stats if row["organelle"] == organelle]
        complete_selected = [row for row in complete_stats if row["organelle"] == organelle]
        labels = [f"{row['component']}:{row['technical_variable'][:4]}" for row in selected]
        values = [float(row["spearman_rho"]) for row in selected]
        complete_values = [float(row["spearman_rho"]) for row in complete_selected]
        x = np.arange(len(values))
        axis.bar(x - 0.18, values, width=0.36, color="#4c78a8", label="Canonical PCA")
        axis.bar(x + 0.18, complete_values, width=0.36, color="#f58518", label="Fully called-site PCA")
        axis.axhline(0, color="black", linewidth=0.8)
        summary = complete_summary[organelle]
        axis.set(
            xticks=x,
            xticklabels=labels,
            title=(
                f"{organelle.capitalize()} PC–QC tests\nComplete-site Procrustes r={float(summary['protest_r']):.3f} ({summary['status']})"
            ),
            ylabel="Spearman ρ",
        )
        axis.tick_params(axis="x", labelrotation=55, labelsize=7)
        axis.legend(fontsize=7)
    fig.suptitle(
        "Supplementary Figure 3 — Genotypes and technical confounders\n"
        "Associations do not distinguish biological divergence from reference-mapping bias; residual NUMT/NUPT ambiguity remains"
    )
    fig.tight_layout()
    return _save(fig, prefix)


def _organelle_comparison(root: Path, run_id: str, prefix: Path) -> list[Path]:
    base = root / f"supplementary_analysis/results/comparative/{run_id}/organelle_comparison"
    tangle = read_tsv(base / "tanglegram_271_tip_mapping.tsv")
    fst = read_tsv(base / "common_pair_fst_agreement.tsv")
    rf = read_tsv(base / "supported_unrooted_rf.tsv")[0]
    fig, axes = plt.subplots(1, 2, figsize=(12, 7))
    for row in tangle:
        axes[0].plot(
            [0, 1], [int(row["chloroplast_tip_order"]), int(row["mitochondria_tip_order"])], color="#777777", alpha=0.18, linewidth=0.5
        )
    axes[0].set(
        xticks=[0, 1], xticklabels=["chloroplast", "mitochondria"], ylabel="Unrooted display tip order", title="All 271 shared samples"
    )
    cp = [float(row["chloroplast_hudson_fst"]) for row in fst]
    mt = [float(row["mitochondria_hudson_fst"]) for row in fst]
    axes[1].scatter(cp, mt, s=12, alpha=0.55, color="#355f8d")
    low, high = min(cp + mt), max(cp + mt)
    axes[1].plot([low, high], [low, high], linestyle="--", color="black", linewidth=1)
    axes[1].set(
        xlabel="Chloroplast signed Hudson FST",
        ylabel="Mitochondrial signed Hudson FST",
        title=f"Common population pairs (ρ={float(fst[0]['global_spearman_rho']):.2f})",
    )
    fig.text(
        0.5,
        0.01,
        f"Supported unrooted RF = {rf['rf_numerator']}/{rf['rf_denominator']} = "
        f"{float(rf['normalized_unrooted_rf']):.3f}; computed on 229 mt representatives",
        ha="center",
        fontsize=9,
    )
    fig.suptitle("Supplementary Figure 4 — Chloroplast–mitochondrial comparison")
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    return _save(fig, prefix)


def _population_diversity(root: Path, run_id: str, prefix: Path) -> list[Path]:
    base = root / f"supplementary_analysis/results/comparative/{run_id}/population_diversity"
    site = read_tsv(base / "chloroplast_146_site_resampling.tsv")
    pi = read_tsv(base / "population_pi_n4_resampling.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    values = [int(row["cp_multi_population_haplotypes"]) for row in site]
    observed = int(site[0]["observed_mt_multi_population_haplotypes"])
    axes[0].hist(values, bins=np.arange(min(values), max(values) + 2) - 0.5, color="#27864a", alpha=0.8)
    axes[0].axvline(observed, color="#8a4f9e", linewidth=2, label=f"Observed mt = {observed}")
    axes[0].legend()
    axes[0].set(title="146-site marker-count sensitivity", xlabel="Multi-population haplotypes", ylabel="1,000 cp draws")
    grouped: dict[str, list[float]] = {}
    for row in pi:
        grouped.setdefault(row["population"], []).append(float(row["nucleotide_diversity"]))
    populations = sorted(grouped, key=lambda key: np.median(grouped[key]), reverse=True)
    selected = populations[:12]
    axes[1].boxplot([grouped[key] for key in selected], tick_labels=selected, showfliers=False)
    axes[1].tick_params(axis="x", labelrotation=60, labelsize=7)
    axes[1].set(title="Sample-size-standardized chloroplast π", ylabel="π (n=4)")
    fig.suptitle("Supplementary Figure 5 — Population diversity (DUSE excluded from inference)")
    fig.tight_layout()
    return _save(fig, prefix)


def _coordinate(root: Path, run_id: str, prefix: Path) -> list[Path]:
    rows = read_tsv(root / f"supplementary_analysis/results/comparative/{run_id}/genome_coordinate_windows.tsv")
    fig, axes = plt.subplots(2, 1, figsize=(12, 7))
    for axis, organelle in zip(axes, ("chloroplast", "mitochondria"), strict=True):
        selected = [row for row in rows if row["organelle"] == organelle]
        x = [(int(row["start_0based"]) + int(row["end_0based_exclusive"])) / 2000 for row in selected]
        variation = [float(row["variation_per_callable_kb"]) if row["variation_per_callable_kb"] != "nan" else np.nan for row in selected]
        depth = [float(row["mean_filtered_depth"]) for row in selected]
        axis.plot(x, variation, color=COLORS[organelle], marker="o", markersize=3, label="SNPs/callable kb")
        twin = axis.twinx()
        twin.plot(x, depth, color="#777777", alpha=0.55, label="Mean filtered depth")
        axis.set(title=organelle.capitalize(), xlabel="Genome coordinate (kb)", ylabel="Variation per callable kb")
        twin.set_ylabel("Mean filtered depth")
    fig.suptitle("Supplementary Figure 6 — Non-overlapping 5 kb genome-coordinate tracks")
    fig.tight_layout()
    return _save(fig, prefix)


def _presentation_replacements(root: Path, run_id: str, output_dir: Path) -> tuple[list[Path], list[dict[str, str]]]:
    outputs: list[Path] = []
    rows: list[dict[str, str]] = []
    # Compact, explicitly demoted summaries; these do not count as supplementary families.
    for organelle in ("chloroplast", "mitochondria"):
        variance = read_tsv(root / f"canonical_publication/results/pca/publication-20260817/{organelle}.variance.tsv")
        fig, axis = plt.subplots(figsize=(6, 4))
        axis.plot(range(1, len(variance) + 1), [100 * float(row["explained_variance_ratio"]) for row in variance], marker="o")
        axis.set(
            title=f"{organelle.capitalize()} PCA scree (presentation replacement)",
            xlabel="Principal component",
            ylabel="Explained variance (%)",
        )
        paths = _save(fig, output_dir / f"replacement_{organelle}_pca_scree")
        outputs.extend(paths)
        rows.extend(
            {
                "replacement_id": f"pca_scree_{organelle}",
                "format": path.suffix[1:],
                "path": path.relative_to(root).as_posix(),
                "supplement_family": "no",
            }
            for path in paths
        )
    fst = read_tsv(root / f"supplementary_analysis/results/comparative/{run_id}/organelle_comparison/common_pair_fst_agreement.tsv")
    fig, axis = plt.subplots(figsize=(6, 4))
    for organelle, field in (("chloroplast", "chloroplast_hudson_fst"), ("mitochondria", "mitochondria_hudson_fst")):
        values = np.sort([float(row[field]) for row in fst])
        axis.step(values, np.arange(1, len(values) + 1) / len(values), label=organelle, color=COLORS[organelle])
    axis.set(title="Signed Hudson FST ECDF (presentation replacement)", xlabel="FST", ylabel="Cumulative fraction")
    axis.legend()
    paths = _save(fig, output_dir / "replacement_fst_ecdf")
    outputs.extend(paths)
    rows.extend(
        {"replacement_id": "fst_ecdf", "format": path.suffix[1:], "path": path.relative_to(root).as_posix(), "supplement_family": "no"}
        for path in paths
    )
    k_rows = read_tsv(root / "canonical_publication/results/supplement/publication-20260817/admixture/chloroplast/k_summary.tsv")
    fig, axis = plt.subplots(figsize=(6, 4))
    cv_field = next(field for field in k_rows[0] if "cv" in field.lower() and field.lower() != "cv_sd")
    k_field = next(field for field in k_rows[0] if field.lower() == "k")
    axis.plot([int(row[k_field]) for row in k_rows], [float(row[cv_field]) for row in k_rows], marker="o")
    axis.set(title="ADMIXTURE — demoted sensitivity view", xlabel="K", ylabel="Cross-validation error")
    axis.text(0.02, 0.02, "Linked haploid organelle markers; not primary population inference", transform=axis.transAxes, fontsize=8)
    paths = _save(fig, output_dir / "replacement_admixture_demoted")
    outputs.extend(paths)
    rows.extend(
        {
            "replacement_id": "admixture_demoted",
            "format": path.suffix[1:],
            "path": path.relative_to(root).as_posix(),
            "supplement_family": "no",
        }
        for path in paths
    )
    return outputs, rows


def render_all_figures(root: Path, run_id: str) -> list[Path]:
    output_dir = root / f"supplementary_analysis/reports/figures/{run_id}"
    functions = (_robustness, _phylogenetic_information, _technical, _organelle_comparison, _population_diversity, _coordinate)
    outputs: list[Path] = []
    manifest_rows: list[dict[str, str]] = []
    for index, (family, function) in enumerate(zip(FIGURE_FAMILIES, functions, strict=True), 1):
        paths = function(root, run_id, output_dir / f"figure_S{index}_{family}")
        outputs.extend(paths)
        manifest_rows.extend(
            {"figure_id": f"S{index}", "family": family, "format": path.suffix[1:], "path": path.relative_to(root).as_posix()}
            for path in paths
        )
    validate_figure_manifest(manifest_rows)
    manifest = output_dir / "supplementary_figure_manifest.tsv"
    write_tsv(manifest, manifest_rows, ["figure_id", "family", "format", "path"], root)
    replacement_outputs, replacement_rows = _presentation_replacements(root, run_id, output_dir)
    replacement_manifest = output_dir / "presentation_replacement_manifest.tsv"
    write_tsv(replacement_manifest, replacement_rows, ["replacement_id", "format", "path", "supplement_family"], root)
    return [*outputs, manifest, *replacement_outputs, replacement_manifest]
