"""Run cpDNA and mtDNA admixture-style clustering.

This stage converts filtered haploid organelle SNP
alignments into pseudo-diploid homozygous PLINK PED/MAP inputs, runs ADMIXTURE
across a fixed K range with cross-validation, picks the lowest-CV K separately
for cpDNA and mtDNA, and renders structure-style plots.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, stdev

from dudleya_organelle_alignment_pipeline.pilot_alignment import shlex_join
from dudleya_organelle_alignment_pipeline.pca_analysis import (
    choose_plot_group,
    read_fasta,
    read_sample_metadata,
)
from dudleya_organelle_alignment_pipeline.variant_calling import labeled_output_name


DEFAULT_SNP_ALIGNMENT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/10_snp_alignment"
)
DEFAULT_METADATA_PATH = Path(
    "dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/included_samples.tsv"
)
DEFAULT_OUTPUT_DIR = Path("dudleya_organelle_alignment_pipeline/results/16_admixture")
DEFAULT_RUN_LABEL = "primary"
DEFAULT_MIN_K = 1
DEFAULT_MAX_K = 8
DEFAULT_THREADS = 4
DEFAULT_SEED = 20260707
BASES = {"A", "C", "G", "T"}


class AdmixtureAnalysisError(RuntimeError):
    """Raised when this stage cannot run safely."""


@dataclass(frozen=True)
class AdmixtureInput:
    organelle: str
    track_id: str
    sample_count: int
    alignment_sites: int
    missing_bases: int
    alignment_fasta_path: Path
    site_table_path: Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def require_admixture() -> str:
    executable = shutil.which("admixture")
    if executable is None:
        raise AdmixtureAnalysisError(
            "Missing required tool: admixture. Activate the pipeline environment first."
        )
    return executable


def require_plink() -> str:
    executable = shutil.which("plink")
    if executable is None:
        raise AdmixtureAnalysisError(
            "Missing required tool: plink. Activate the pipeline environment first."
        )
    return executable


def read_admixture_inputs(
    snp_alignment_dir: Path = DEFAULT_SNP_ALIGNMENT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
) -> list[AdmixtureInput]:
    summary_path = snp_alignment_dir / labeled_output_name(
        "snp_alignment_summary.tsv",
        run_label,
    )
    inputs: list[AdmixtureInput] = []
    for row in read_tsv(summary_path):
        fasta_path = Path(row["alignment_fasta_path"])
        site_table_path = Path(row["site_table_path"])
        if not fasta_path.exists() or fasta_path.stat().st_size == 0:
            raise AdmixtureAnalysisError(f"Missing SNP alignment FASTA: {fasta_path}")
        inputs.append(
            AdmixtureInput(
                organelle=row["organelle"],
                track_id=row["track_id"],
                sample_count=int(row["sample_count"]),
                alignment_sites=int(row["alignment_sites"]),
                missing_bases=int(row["missing_bases"]),
                alignment_fasta_path=fasta_path,
                site_table_path=site_table_path,
            )
        )
    if not inputs:
        raise AdmixtureAnalysisError(f"No SNP alignment rows found in {summary_path}")
    return inputs


def write_pseudo_diploid_ped_map(
    admixture_input: AdmixtureInput,
    output_dir: Path,
    run_label: str,
) -> tuple[Path, Path, list[str], list[str]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    records = read_fasta(admixture_input.alignment_fasta_path)
    prefix = output_dir / f"{admixture_input.organelle}.{run_label}.pseudo_diploid"
    ped_path = Path(f"{prefix}.ped")
    map_path = Path(f"{prefix}.map")
    excluded_path = output_dir / (
        f"{admixture_input.organelle}.{run_label}.pseudo_diploid.excluded_samples.tsv"
    )
    included_sample_ids: list[str] = []
    excluded_sample_ids: list[str] = []

    with map_path.open("w") as map_handle:
        for site_index in range(admixture_input.alignment_sites):
            marker_id = f"{admixture_input.organelle}_snp_{site_index + 1}"
            map_handle.write(f"1\t{marker_id}\t0\t{site_index + 1}\n")

    with ped_path.open("w") as ped_handle:
        for sample_id, sequence in records:
            if all(base not in BASES for base in sequence):
                excluded_sample_ids.append(sample_id)
                continue
            included_sample_ids.append(sample_id)
            genotype_fields: list[str] = []
            for base in sequence:
                allele = base if base in BASES else "0"
                genotype_fields.extend([allele, allele])
            ped_handle.write(
                " ".join([sample_id, sample_id, "0", "0", "0", "-9", *genotype_fields])
                + "\n"
            )
    with excluded_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample_id", "organelle", "reason"],
            delimiter="\t",
        )
        writer.writeheader()
        for sample_id in excluded_sample_ids:
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "organelle": admixture_input.organelle,
                    "reason": "all_snp_genotypes_missing",
                }
            )
    if not included_sample_ids:
        raise AdmixtureAnalysisError(
            f"No informative {admixture_input.organelle} samples remain for ADMIXTURE"
        )
    return ped_path, map_path, included_sample_ids, excluded_sample_ids


def build_admixture_command(
    admixture_executable: str,
    genotype_path: Path,
    k: int,
    seed: int,
    threads: int,
) -> list[str]:
    return [
        admixture_executable,
        "--cv",
        f"--seed={seed}",
        f"-j{threads}",
        genotype_path.name,
        str(k),
    ]


def build_plink_make_bed_command(plink_executable: str, ped_prefix: Path) -> list[str]:
    return [
        plink_executable,
        "--file",
        ped_prefix.name,
        "--make-bed",
        "--out",
        ped_prefix.name,
    ]


def run_plink_make_bed(
    plink_executable: str,
    ped_path: Path,
    output_dir: Path,
    force: bool,
) -> tuple[Path, str]:
    ped_prefix = ped_path.with_suffix("")
    bed_path = Path(f"{ped_prefix}.bed")
    bim_path = Path(f"{ped_prefix}.bim")
    fam_path = Path(f"{ped_prefix}.fam")
    command = build_plink_make_bed_command(plink_executable, ped_prefix)
    if bed_path.exists() and bim_path.exists() and fam_path.exists() and not force:
        return bed_path, "outputs already present; pass --force to regenerate"
    completed = subprocess.run(
        command,
        cwd=output_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path = output_dir / f"{ped_prefix.name}.plink_make_bed.log"
    log_path.write_text(completed.stdout)
    if completed.returncode:
        raise AdmixtureAnalysisError(f"PLINK conversion failed; see {log_path}")
    return bed_path, shlex_join(command)


def parse_cv_error(log_text: str) -> float:
    match = re.search(r"CV error \(K=\d+\):\s*([0-9.eE+-]+)", log_text)
    if match is None:
        raise AdmixtureAnalysisError("Could not parse ADMIXTURE CV error from log")
    return float(match.group(1))


def summarize_replicate_stability(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        key = (row["organelle"], row["k"])
        grouped.setdefault(key, []).append(float(row["cv_error"]))
    best_by_organelle: dict[str, str] = {}
    for organelle in sorted({organelle for organelle, _ in grouped}):
        organelle_items = [
            (k, mean(values)) for (group_organelle, k), values in grouped.items()
            if group_organelle == organelle
        ]
        best_by_organelle[organelle] = min(organelle_items, key=lambda item: item[1])[0]
    summary: list[dict[str, str]] = []
    for organelle, k in sorted(grouped, key=lambda item: (item[0], int(item[1]))):
        values = grouped[(organelle, k)]
        summary.append(
            {
                "organelle": organelle,
                "k": k,
                "replicate_count": str(len(values)),
                "mean_cv_error": f"{mean(values):.8f}",
                "sd_cv_error": f"{(stdev(values) if len(values) > 1 else 0.0):.8f}",
                "min_cv_error": f"{min(values):.8f}",
                "max_cv_error": f"{max(values):.8f}",
                "is_best_mean_k": "yes" if k == best_by_organelle[organelle] else "no",
            }
        )
    return summary


def run_admixture_for_k(
    admixture_executable: str,
    genotype_path: Path,
    output_dir: Path,
    organelle: str,
    k: int,
    replicate: int,
    seed: int,
    threads: int,
    force: bool,
) -> dict[str, str]:
    replicate_label = f".rep{replicate}"
    base_q_path = output_dir / f"{genotype_path.stem}.{k}.Q"
    base_p_path = output_dir / f"{genotype_path.stem}.{k}.P"
    q_path = output_dir / f"{genotype_path.stem}.{k}{replicate_label}.Q"
    p_path = output_dir / f"{genotype_path.stem}.{k}{replicate_label}.P"
    log_path = output_dir / f"{organelle}.K{k}.rep{replicate}.admixture.log"
    command = build_admixture_command(
        admixture_executable=admixture_executable,
        genotype_path=genotype_path,
        k=k,
        seed=seed + replicate - 1,
        threads=threads,
    )
    if q_path.exists() and p_path.exists() and log_path.exists() and not force:
        log_text = log_path.read_text(errors="replace")
    else:
        completed = subprocess.run(
            command,
            cwd=output_dir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        log_text = completed.stdout
        log_path.write_text(log_text)
        if completed.returncode:
            raise AdmixtureAnalysisError(
                f"ADMIXTURE failed for {organelle} K={k}; see {log_path}"
            )
        if not base_q_path.exists() or not base_p_path.exists():
            raise AdmixtureAnalysisError(
                f"ADMIXTURE did not write expected Q/P files for {organelle} K={k}"
            )
        base_q_path.replace(q_path)
        base_p_path.replace(p_path)
    cv_error = parse_cv_error(log_text)
    return {
        "organelle": organelle,
        "k": str(k),
        "replicate": str(replicate),
        "cv_error": f"{cv_error:.8f}",
        "q_path": q_path.as_posix(),
        "p_path": p_path.as_posix(),
        "log_path": log_path.as_posix(),
        "command": shlex_join(command),
    }


def read_q_matrix(q_path: Path) -> list[list[float]]:
    matrix: list[list[float]] = []
    with q_path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                matrix.append([float(value) for value in stripped.split()])
    return matrix


def write_q_table_and_plot(
    admixture_input: AdmixtureInput,
    q_path: Path,
    metadata_path: Path,
    output_dir: Path,
    run_label: str,
    best_k: int,
    sample_ids: list[str] | None = None,
) -> tuple[Path, Path, Path, Path]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise AdmixtureAnalysisError("Missing matplotlib for ADMIXTURE plots") from exc

    if sample_ids is None:
        sample_ids = [
            sample_id for sample_id, _ in read_fasta(admixture_input.alignment_fasta_path)
        ]
    q_matrix = read_q_matrix(q_path)
    if len(q_matrix) != len(sample_ids):
        raise AdmixtureAnalysisError(f"Q row count does not match samples for {q_path}")
    metadata = read_sample_metadata(metadata_path)
    ordered = sorted(
        enumerate(sample_ids),
        key=lambda item: (
            choose_plot_group(metadata.get(item[1], {})),
            item[1],
        ),
    )

    table_path = output_dir / f"{admixture_input.organelle}.{run_label}.bestK{best_k}.q.tsv"
    q_fields = [f"cluster_{index + 1}" for index in range(best_k)]
    rows: list[dict[str, str]] = []
    for original_index, sample_id in ordered:
        meta = metadata.get(sample_id, {})
        row = {
            "sample_id": sample_id,
            "organelle": admixture_input.organelle,
            "best_k": str(best_k),
            "species": meta.get("species", ""),
            "popcode": meta.get("popcode", ""),
            "population_name": meta.get("population_name", ""),
            "plot_group": choose_plot_group(meta),
        }
        for q_index, q_value in enumerate(q_matrix[original_index]):
            row[f"cluster_{q_index + 1}"] = f"{q_value:.8f}"
        rows.append(row)
    write_tsv(
        table_path,
        rows,
        [
            "sample_id",
            "organelle",
            "best_k",
            "species",
            "popcode",
            "population_name",
            "plot_group",
            *q_fields,
        ],
    )

    fig_width = 14
    fig_height = max(4, min(8, len(sample_ids) * 0.018 + 3))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    bottoms = [0.0] * len(rows)
    colors = plt.get_cmap("tab20").colors
    x_values = list(range(len(rows)))
    for cluster_index in range(best_k):
        values = [float(row[f"cluster_{cluster_index + 1}"]) for row in rows]
        ax.bar(
            x_values,
            values,
            bottom=bottoms,
            width=1.0,
            color=colors[cluster_index % len(colors)],
            edgecolor="none",
            label=f"K{cluster_index + 1}",
        )
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
    ax.set_title(f"{admixture_input.organelle} ADMIXTURE-style clustering, K={best_k}")
    ax.set_ylabel("Assignment proportion")
    ax.set_xlabel("Samples sorted by metadata group")
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.legend(ncol=min(best_k, 8), loc="upper center", bbox_to_anchor=(0.5, -0.12))

    prefix = output_dir / f"{admixture_input.organelle}.{run_label}.bestK{best_k}.structure"
    png_path = Path(f"{prefix}.png")
    pdf_path = Path(f"{prefix}.pdf")
    svg_path = Path(f"{prefix}.svg")
    fig.savefig(png_path, dpi=200)
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)
    return table_path, png_path, pdf_path, svg_path


def write_cv_plot(
    output_dir: Path,
    organelle: str,
    run_label: str,
    rows: list[dict[str, str]],
) -> Path:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise AdmixtureAnalysisError("Missing matplotlib for ADMIXTURE CV plot") from exc

    summary_rows = summarize_replicate_stability(rows)
    ordered = sorted(
        [row for row in summary_rows if row["organelle"] == organelle],
        key=lambda row: int(row["k"]),
    )
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.plot(
        [int(row["k"]) for row in ordered],
        [float(row["mean_cv_error"]) for row in ordered],
        marker="o",
        linewidth=1.5,
    )
    ax.set_title(f"{organelle} ADMIXTURE CV error")
    ax.set_xlabel("K")
    ax.set_ylabel("Mean cross-validation error")
    ax.grid(alpha=0.25, linewidth=0.5)
    png_path = output_dir / f"{organelle}.{run_label}.admixture_cv.png"
    fig.savefig(png_path, dpi=200)
    plt.close(fig)
    return png_path


def write_admixture_outputs(
    output_dir: Path,
    rows: list[dict[str, str]],
    run_label: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / labeled_output_name("admixture_summary.tsv", run_label),
        rows,
        [
            "organelle",
            "track_id",
            "k",
            "replicate",
            "cv_error",
            "is_best_k",
            "is_best_mean_k",
            "mean_cv_error",
            "sd_cv_error",
            "replicate_count",
            "excluded_sample_count",
            "q_path",
            "p_path",
            "log_path",
            "best_q_table_path",
            "structure_png_path",
            "structure_pdf_path",
            "structure_svg_path",
            "cv_plot_path",
            "plink_command",
            "command",
        ],
    )
    write_admixture_report(
        output_dir / labeled_output_name("admixture_report.md", run_label),
        rows,
        run_label,
    )


def write_admixture_report(path: Path, rows: list[dict[str, str]], run_label: str) -> None:
    label = run_label or "full"
    lines = [
        "# Admixture-Style Clustering",
        "",
        "This step runs ADMIXTURE on cpDNA and mtDNA SNP alignments separately.",
        "Because ADMIXTURE is a diploid-oriented tool, haploid organelle calls",
        "are encoded as pseudo-diploid homozygous genotypes. These plots should",
        "be interpreted as organelle haplotype clustering, not nuclear admixture.",
        "",
        "## Run",
        "",
        f"- Run label: `{label}`",
        "- K selection: lowest mean cross-validation error among tested K values",
        "- Haploid encoding: each called base is duplicated; missing calls are `0 0`",
        "",
        "## Results",
        "",
    ]
    for organelle in sorted({row["organelle"] for row in rows}):
        organelle_rows = sorted(
            [row for row in rows if row["organelle"] == organelle],
            key=lambda row: int(row["k"]),
        )
        best = next(row for row in organelle_rows if row["is_best_k"] == "yes")
        lines.extend(
            [
                f"### {organelle}",
                "",
                f"- Track: `{best['track_id']}`",
                f"- Best K: {best['k']}",
                f"- Mean CV error at best K: {best['mean_cv_error']}",
                f"- CV-error SD at best K: {best['sd_cv_error']}",
                f"- Replicates per K: {best['replicate_count']}",
                f"- Structure plot: `{best['structure_png_path']}`",
                f"- Best-K Q table: `{best['best_q_table_path']}`",
                f"- CV plot: `{best['cv_plot_path']}`",
                "",
                "| K | Mean CV error | SD | Replicates | Best |",
                "|---|---:|---:|---:|---|",
            ]
        )
        seen_k = set()
        for row in organelle_rows:
            if row["k"] in seen_k:
                continue
            seen_k.add(row["k"])
            lines.append(
                f"| {row['k']} | {row['mean_cv_error']} | {row['sd_cv_error']} | "
                f"{row['replicate_count']} | {row['is_best_mean_k']} |"
            )
        lines.append("")
    path.write_text("\n".join(lines))


def run_admixture_analysis(
    snp_alignment_dir: Path = DEFAULT_SNP_ALIGNMENT_DIR,
    metadata_path: Path = DEFAULT_METADATA_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
    min_k: int = DEFAULT_MIN_K,
    max_k: int = DEFAULT_MAX_K,
    threads: int = DEFAULT_THREADS,
    seed: int = DEFAULT_SEED,
    replicates: int = 1,
    force: bool = False,
) -> list[dict[str, str]]:
    if min_k < 1 or max_k < min_k:
        raise AdmixtureAnalysisError("Invalid K range")
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/dudleya_matplotlib")
    output_dir.mkdir(parents=True, exist_ok=True)
    admixture_executable = require_admixture()
    plink_executable = require_plink()
    inputs = read_admixture_inputs(snp_alignment_dir=snp_alignment_dir, run_label=run_label)
    all_rows: list[dict[str, str]] = []
    for admixture_input in inputs:
        ped_path, _, included_sample_ids, excluded_sample_ids = write_pseudo_diploid_ped_map(
            admixture_input,
            output_dir,
            run_label,
        )
        bed_path, plink_command = run_plink_make_bed(
            plink_executable=plink_executable,
            ped_path=ped_path,
            output_dir=output_dir,
            force=force,
        )
        k_max_for_input = min(
            max_k,
            len(included_sample_ids) - 1,
            admixture_input.alignment_sites,
        )
        if k_max_for_input < min_k:
            raise AdmixtureAnalysisError(
                f"Not enough informative {admixture_input.organelle} samples for K range"
            )
        k_rows: list[dict[str, str]] = []
        for k in range(min_k, k_max_for_input + 1):
            for replicate in range(1, replicates + 1):
                row = run_admixture_for_k(
                    admixture_executable=admixture_executable,
                    genotype_path=bed_path,
                    output_dir=output_dir,
                    organelle=admixture_input.organelle,
                    k=k,
                    replicate=replicate,
                    seed=seed,
                    threads=threads,
                    force=force,
                )
                row["track_id"] = admixture_input.track_id
                row["plink_command"] = plink_command
                row["excluded_sample_count"] = str(len(excluded_sample_ids))
                k_rows.append(row)
        stability_rows = summarize_replicate_stability(k_rows)
        stability_by_k = {row["k"]: row for row in stability_rows if row["organelle"] == admixture_input.organelle}
        best_k_from_mean = next(row["k"] for row in stability_rows if row["organelle"] == admixture_input.organelle and row["is_best_mean_k"] == "yes")
        best_row = min(
            [row for row in k_rows if row["k"] == best_k_from_mean],
            key=lambda row: float(row["cv_error"]),
        )
        best_k = int(best_row["k"])
        q_table_path, structure_png, structure_pdf, structure_svg = write_q_table_and_plot(
            admixture_input=admixture_input,
            q_path=Path(best_row["q_path"]),
            metadata_path=metadata_path,
            output_dir=output_dir,
            run_label=run_label,
            best_k=best_k,
            sample_ids=included_sample_ids,
        )
        cv_plot_path = write_cv_plot(output_dir, admixture_input.organelle, run_label, k_rows)
        for row in k_rows:
            stability = stability_by_k[row["k"]]
            row["mean_cv_error"] = stability["mean_cv_error"]
            row["sd_cv_error"] = stability["sd_cv_error"]
            row["replicate_count"] = stability["replicate_count"]
            row["is_best_mean_k"] = stability["is_best_mean_k"]
            row["is_best_k"] = "yes" if row["k"] == str(best_k) and row["replicate"] == best_row["replicate"] else "no"
            row["best_q_table_path"] = q_table_path.as_posix() if row["is_best_k"] == "yes" else ""
            row["structure_png_path"] = structure_png.as_posix() if row["is_best_k"] == "yes" else ""
            row["structure_pdf_path"] = structure_pdf.as_posix() if row["is_best_k"] == "yes" else ""
            row["structure_svg_path"] = structure_svg.as_posix() if row["is_best_k"] == "yes" else ""
            row["cv_plot_path"] = cv_plot_path.as_posix()
        all_rows.extend(k_rows)
    write_admixture_outputs(output_dir, all_rows, run_label=run_label)
    return all_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run cpDNA/mtDNA ADMIXTURE-style clustering."
    )
    parser.add_argument("--snp-alignment-dir", type=Path, default=DEFAULT_SNP_ALIGNMENT_DIR)
    parser.add_argument("--metadata-path", type=Path, default=DEFAULT_METADATA_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    parser.add_argument("--min-k", type=int, default=DEFAULT_MIN_K)
    parser.add_argument("--max-k", type=int, default=DEFAULT_MAX_K)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--replicates", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    rows = run_admixture_analysis(
        snp_alignment_dir=args.snp_alignment_dir,
        metadata_path=args.metadata_path,
        output_dir=args.output_dir,
        run_label=args.run_label,
        min_k=args.min_k,
        max_k=args.max_k,
        threads=args.threads,
        seed=args.seed,
        replicates=args.replicates,
        force=args.force,
    )
    for organelle in sorted({row["organelle"] for row in rows}):
        best = next(
            row
            for row in rows
            if row["organelle"] == organelle and row["is_best_k"] == "yes"
        )
        print(
            f"{organelle}: best K={best['k']} with mean CV error {best['mean_cv_error']} "
            f"at {best['structure_png_path']}"
        )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
