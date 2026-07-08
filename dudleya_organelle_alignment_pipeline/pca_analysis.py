"""Build cpDNA and mtDNA PCA visualizations from haploid SNP alignments.

This stage consumes the SNP-only FASTA
alignments and writes PCA coordinates, variance summaries, static figures, and
a short report for cpDNA and mtDNA separately.
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from dudleya_organelle_alignment_pipeline.variant_calling import labeled_output_name


DEFAULT_SNP_ALIGNMENT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/10_snp_alignment"
)
DEFAULT_METADATA_PATH = Path(
    "dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/included_samples.tsv"
)
DEFAULT_OUTPUT_DIR = Path("dudleya_organelle_alignment_pipeline/results/15_pca")
DEFAULT_RUN_LABEL = "primary"
BASES = ("A", "C", "G", "T")


class PcaAnalysisError(RuntimeError):
    """Raised when PCA inputs cannot be converted into a figure-ready matrix."""


@dataclass(frozen=True)
class PcaInput:
    organelle: str
    track_id: str
    sample_count: int
    alignment_sites: int
    missing_bases: int
    alignment_fasta_path: Path
    site_table_path: Path


@dataclass(frozen=True)
class PcaResult:
    organelle: str
    track_id: str
    sample_count: int
    alignment_sites: int
    retained_sites: int
    missing_bases: int
    pc1_variance: float
    pc2_variance: float
    alignment_fasta_path: Path
    coordinates_path: Path
    variance_path: Path
    png_path: Path
    pdf_path: Path
    svg_path: Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_pca_inputs(
    snp_alignment_dir: Path = DEFAULT_SNP_ALIGNMENT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
) -> list[PcaInput]:
    summary_path = snp_alignment_dir / labeled_output_name(
        "snp_alignment_summary.tsv",
        run_label,
    )
    inputs: list[PcaInput] = []
    for row in read_tsv(summary_path):
        fasta_path = Path(row["alignment_fasta_path"])
        site_table_path = Path(row["site_table_path"])
        if not fasta_path.exists() or fasta_path.stat().st_size == 0:
            raise PcaAnalysisError(f"Missing SNP alignment FASTA: {fasta_path}")
        inputs.append(
            PcaInput(
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
        raise PcaAnalysisError(f"No SNP alignment rows found in {summary_path}")
    return inputs


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    current_name: str | None = None
    sequence_parts: list[str] = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records.append((current_name, "".join(sequence_parts).upper()))
                current_name = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)
    if current_name is not None:
        records.append((current_name, "".join(sequence_parts).upper()))
    if not records:
        raise PcaAnalysisError(f"No FASTA records found in {path}")
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise PcaAnalysisError(f"FASTA records have inconsistent lengths in {path}")
    return records


def build_haploid_snp_matrix(fasta_path: Path) -> tuple[np.ndarray, list[str], int]:
    records = read_fasta(fasta_path)
    sample_ids = [sample_id for sample_id, _ in records]
    sequences = [sequence for _, sequence in records]
    site_count = len(sequences[0])
    columns: list[np.ndarray] = []

    for site_index in range(site_count):
        site_bases = [sequence[site_index] for sequence in sequences]
        alleles = sorted({base for base in site_bases if base in BASES})
        if len(alleles) < 2:
            continue
        allele_codes = {base: float(index) for index, base in enumerate(alleles)}
        encoded = np.array(
            [allele_codes.get(base, np.nan) for base in site_bases],
            dtype=float,
        )
        if np.all(np.isnan(encoded)):
            continue
        site_mean = float(np.nanmean(encoded))
        encoded[np.isnan(encoded)] = site_mean
        columns.append(encoded)

    if not columns:
        raise PcaAnalysisError(f"No polymorphic SNP columns retained from {fasta_path}")
    matrix = np.column_stack(columns)
    return matrix, sample_ids, len(columns)


def read_sample_metadata(metadata_path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(metadata_path)
    metadata: dict[str, dict[str, str]] = {}
    for row in rows:
        metadata[row["sample_id"]] = row
    return metadata


def choose_plot_group(row: dict[str, str]) -> str:
    species = row.get("species", "").strip()
    popcode = row.get("popcode", "").strip()
    naming_profile = row.get("naming_profile", "").strip()
    if species and popcode:
        return f"{species}_{popcode}"
    if popcode:
        return popcode
    if species:
        return species
    return naming_profile or "unresolved"


def run_pca(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    try:
        from sklearn.decomposition import PCA
    except ImportError as exc:
        raise PcaAnalysisError(
            "Missing scikit-learn. Activate the pipeline environment or rerun the tool audit."
        ) from exc

    n_components = min(2, matrix.shape[0], matrix.shape[1])
    if n_components < 2:
        raise PcaAnalysisError("PCA requires at least two samples and two retained SNP sites")
    centered = matrix - matrix.mean(axis=0)
    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(centered)
    return coordinates, pca.explained_variance_ratio_


def write_pca_tables(
    output_dir: Path,
    pca_input: PcaInput,
    sample_ids: list[str],
    coordinates: np.ndarray,
    variance: np.ndarray,
    retained_sites: int,
    metadata: dict[str, dict[str, str]],
    run_label: str,
) -> tuple[Path, Path]:
    prefix = output_dir / f"{pca_input.organelle}.{run_label}.pca"
    coordinates_path = Path(f"{prefix}.coordinates.tsv")
    variance_path = Path(f"{prefix}.variance.tsv")

    coordinate_rows: list[dict[str, str]] = []
    for index, sample_id in enumerate(sample_ids):
        row = metadata.get(sample_id, {})
        coordinate_rows.append(
            {
                "sample_id": sample_id,
                "organelle": pca_input.organelle,
                "pc1": f"{coordinates[index, 0]:.8f}",
                "pc2": f"{coordinates[index, 1]:.8f}",
                "species": row.get("species", ""),
                "popcode": row.get("popcode", ""),
                "population_name": row.get("population_name", ""),
                "naming_profile": row.get("naming_profile", ""),
                "plot_group": choose_plot_group(row),
            }
        )
    write_tsv(
        coordinates_path,
        coordinate_rows,
        [
            "sample_id",
            "organelle",
            "pc1",
            "pc2",
            "species",
            "popcode",
            "population_name",
            "naming_profile",
            "plot_group",
        ],
    )

    write_tsv(
        variance_path,
        [
            {
                "organelle": pca_input.organelle,
                "component": "PC1",
                "explained_variance_ratio": f"{variance[0]:.8f}",
                "retained_sites": str(retained_sites),
            },
            {
                "organelle": pca_input.organelle,
                "component": "PC2",
                "explained_variance_ratio": f"{variance[1]:.8f}",
                "retained_sites": str(retained_sites),
            },
        ],
        ["organelle", "component", "explained_variance_ratio", "retained_sites"],
    )
    return coordinates_path, variance_path


def write_pca_plot(
    coordinates_path: Path,
    variance: np.ndarray,
    organelle: str,
    output_dir: Path,
    run_label: str,
) -> tuple[Path, Path, Path]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise PcaAnalysisError(
            "Missing matplotlib. Activate the pipeline environment or rerun the tool audit."
        ) from exc

    rows = read_tsv(coordinates_path)
    groups = sorted({row["plot_group"] for row in rows})
    color_map = plt.get_cmap("tab20")
    colors = {group: color_map(index % 20) for index, group in enumerate(groups)}

    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
    for group in groups:
        group_rows = [row for row in rows if row["plot_group"] == group]
        ax.scatter(
            [float(row["pc1"]) for row in group_rows],
            [float(row["pc2"]) for row in group_rows],
            s=32,
            alpha=0.82,
            label=group,
            color=colors[group],
            edgecolors="none",
        )
    ax.set_title(f"{organelle} PCA")
    ax.set_xlabel(f"PC1 ({variance[0] * 100:.2f}% variance)")
    ax.set_ylabel(f"PC2 ({variance[1] * 100:.2f}% variance)")
    ax.axhline(0, color="#666666", linewidth=0.6, alpha=0.4)
    ax.axvline(0, color="#666666", linewidth=0.6, alpha=0.4)
    ax.grid(alpha=0.18, linewidth=0.5)
    if len(groups) <= 24:
        ax.legend(
            title="Group",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            fontsize=8,
        )

    prefix = output_dir / f"{organelle}.{run_label}.pca"
    png_path = Path(f"{prefix}.png")
    pdf_path = Path(f"{prefix}.pdf")
    svg_path = Path(f"{prefix}.svg")
    fig.savefig(png_path, dpi=200)
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)
    return png_path, pdf_path, svg_path


def run_one_pca(
    pca_input: PcaInput,
    metadata_path: Path,
    output_dir: Path,
    run_label: str = DEFAULT_RUN_LABEL,
) -> PcaResult:
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = read_sample_metadata(metadata_path)
    matrix, sample_ids, retained_sites = build_haploid_snp_matrix(
        pca_input.alignment_fasta_path
    )
    coordinates, variance = run_pca(matrix)
    coordinates_path, variance_path = write_pca_tables(
        output_dir=output_dir,
        pca_input=pca_input,
        sample_ids=sample_ids,
        coordinates=coordinates,
        variance=variance,
        retained_sites=retained_sites,
        metadata=metadata,
        run_label=run_label,
    )
    png_path, pdf_path, svg_path = write_pca_plot(
        coordinates_path=coordinates_path,
        variance=variance,
        organelle=pca_input.organelle,
        output_dir=output_dir,
        run_label=run_label,
    )
    return PcaResult(
        organelle=pca_input.organelle,
        track_id=pca_input.track_id,
        sample_count=len(sample_ids),
        alignment_sites=pca_input.alignment_sites,
        retained_sites=retained_sites,
        missing_bases=pca_input.missing_bases,
        pc1_variance=float(variance[0]),
        pc2_variance=float(variance[1]),
        alignment_fasta_path=pca_input.alignment_fasta_path,
        coordinates_path=coordinates_path,
        variance_path=variance_path,
        png_path=png_path,
        pdf_path=pdf_path,
        svg_path=svg_path,
    )


def write_pca_outputs(output_dir: Path, results: list[PcaResult], run_label: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / labeled_output_name("pca_summary.tsv", run_label),
        [
            {
                "organelle": result.organelle,
                "track_id": result.track_id,
                "sample_count": str(result.sample_count),
                "alignment_sites": str(result.alignment_sites),
                "retained_sites": str(result.retained_sites),
                "missing_bases": str(result.missing_bases),
                "pc1_variance": f"{result.pc1_variance:.8f}",
                "pc2_variance": f"{result.pc2_variance:.8f}",
                "alignment_fasta_path": result.alignment_fasta_path.as_posix(),
                "coordinates_path": result.coordinates_path.as_posix(),
                "variance_path": result.variance_path.as_posix(),
                "png_path": result.png_path.as_posix(),
                "pdf_path": result.pdf_path.as_posix(),
                "svg_path": result.svg_path.as_posix(),
            }
            for result in results
        ],
        [
            "organelle",
            "track_id",
            "sample_count",
            "alignment_sites",
            "retained_sites",
            "missing_bases",
            "pc1_variance",
            "pc2_variance",
            "alignment_fasta_path",
            "coordinates_path",
            "variance_path",
            "png_path",
            "pdf_path",
            "svg_path",
        ],
    )
    write_pca_report(
        output_dir / labeled_output_name("pca_report.md", run_label),
        results=results,
        run_label=run_label,
    )


def write_pca_report(path: Path, results: list[PcaResult], run_label: str) -> None:
    label = run_label or "full"
    lines = [
        "# PCA Visualization",
        "",
        "This step computes cpDNA and mtDNA PCA from the filtered haploid",
        "SNP-only alignments. Missing SNP states are mean-imputed per retained",
        "site before PCA, and plots are colored by available species/population",
        "metadata.",
        "",
        "## Run",
        "",
        f"- Run label: `{label}`",
        "- Input: SNP-only FASTA alignments",
        "- Output formats: coordinates TSV, variance TSV, PNG, PDF, SVG",
        "",
        "## Results",
        "",
    ]
    for result in results:
        lines.extend(
            [
                f"### {result.organelle}",
                "",
                f"- Track: `{result.track_id}`",
                f"- Samples: {result.sample_count}",
                f"- SNP alignment sites: {result.alignment_sites}",
                f"- Retained polymorphic sites: {result.retained_sites}",
                f"- PC1 variance: {result.pc1_variance * 100:.2f}%",
                f"- PC2 variance: {result.pc2_variance * 100:.2f}%",
                f"- Coordinates: `{result.coordinates_path}`",
                f"- PNG: `{result.png_path}`",
                f"- PDF: `{result.pdf_path}`",
                f"- SVG: `{result.svg_path}`",
                "",
            ]
        )
    path.write_text("\n".join(lines))


def run_pca_analysis(
    snp_alignment_dir: Path = DEFAULT_SNP_ALIGNMENT_DIR,
    metadata_path: Path = DEFAULT_METADATA_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
) -> list[PcaResult]:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/dudleya_matplotlib")
    inputs = read_pca_inputs(snp_alignment_dir=snp_alignment_dir, run_label=run_label)
    results = [
        run_one_pca(pca_input, metadata_path, output_dir, run_label=run_label)
        for pca_input in inputs
    ]
    write_pca_outputs(output_dir, results, run_label=run_label)
    return results


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build cpDNA/mtDNA PCA visualizations."
    )
    parser.add_argument("--snp-alignment-dir", type=Path, default=DEFAULT_SNP_ALIGNMENT_DIR)
    parser.add_argument("--metadata-path", type=Path, default=DEFAULT_METADATA_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    results = run_pca_analysis(
        snp_alignment_dir=args.snp_alignment_dir,
        metadata_path=args.metadata_path,
        output_dir=args.output_dir,
        run_label=args.run_label,
    )
    for result in results:
        print(
            f"{result.organelle}: PCA for {result.sample_count} samples and "
            f"{result.retained_sites} retained SNP sites at {result.png_path}"
        )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
