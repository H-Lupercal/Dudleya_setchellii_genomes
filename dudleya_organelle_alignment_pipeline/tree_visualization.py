"""Render cpDNA and mtDNA phylogenetic tree figures.

This is Step 14 of the pipeline. It consumes the Step 11 IQ-TREE Newick
treefiles and writes static PNG, PDF, and SVG figures for review and reports.
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path

from dudleya_organelle_alignment_pipeline.variant_calling import labeled_output_name


DEFAULT_TREE_DIR = Path("dudleya_organelle_alignment_pipeline/results/12_phylogenetic_tree")
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/14_tree_visualization"
)
DEFAULT_RUN_LABEL = "primary"


class TreeVisualizationError(RuntimeError):
    """Raised when tree visualizations cannot be rendered safely."""


@dataclass(frozen=True)
class TreeFigureInput:
    organelle: str
    track_id: str
    sample_count: int
    alignment_sites: int
    model: str
    method: str
    treefile_path: Path


@dataclass(frozen=True)
class TreeFigureResult:
    organelle: str
    track_id: str
    sample_count: int
    alignment_sites: int
    model: str
    method: str
    treefile_path: Path
    png_path: Path
    pdf_path: Path
    svg_path: Path
    tip_count: int
    figure_width: float
    figure_height: float


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tree_figure_inputs(
    tree_dir: Path = DEFAULT_TREE_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
) -> list[TreeFigureInput]:
    summary_path = tree_dir / labeled_output_name(
        "phylogenetic_tree_summary.tsv",
        run_label,
    )
    rows = read_tsv(summary_path)
    inputs: list[TreeFigureInput] = []
    for row in rows:
        treefile_path = Path(row["treefile_path"])
        if not treefile_path.exists() or treefile_path.stat().st_size == 0:
            raise TreeVisualizationError(f"Missing treefile: {treefile_path}")
        inputs.append(
            TreeFigureInput(
                organelle=row["organelle"],
                track_id=row["track_id"],
                sample_count=int(row["sample_count"]),
                alignment_sites=int(row["alignment_sites"]),
                model=row["model"],
                method=row["method"],
                treefile_path=treefile_path,
            )
        )
    if not inputs:
        raise TreeVisualizationError(f"No tree rows found in {summary_path}")
    return inputs


def compute_tree_figure_size(sample_count: int) -> tuple[float, float]:
    width = 14.0
    height = max(6.0, min(80.0, sample_count * 0.16 + 2.0))
    return width, height


def render_tree_figure(
    figure_input: TreeFigureInput,
    output_dir: Path,
    run_label: str = DEFAULT_RUN_LABEL,
) -> TreeFigureResult:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from Bio import Phylo
    except ImportError as exc:
        raise TreeVisualizationError(
            "Missing plotting dependency. Run the tool audit and activate the pipeline environment."
        ) from exc

    output_dir.mkdir(parents=True, exist_ok=True)
    tree = Phylo.read(figure_input.treefile_path, "newick")
    tip_count = len(tree.get_terminals())
    width, height = compute_tree_figure_size(max(figure_input.sample_count, tip_count))
    prefix = output_dir / f"{figure_input.organelle}.{run_label}.iqtree_ml_tree"
    png_path = Path(f"{prefix}.png")
    pdf_path = Path(f"{prefix}.pdf")
    svg_path = Path(f"{prefix}.svg")

    fig = plt.figure(figsize=(width, height), constrained_layout=True)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(
        f"{figure_input.organelle} phylogenetic tree",
        fontsize=14,
        fontweight="bold",
    )
    Phylo.draw(tree, axes=ax, do_show=False)
    ax.set_xlabel("Substitutions per site")
    ax.grid(axis="x", alpha=0.25, linewidth=0.5)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    fig.savefig(png_path, dpi=200)
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)

    return TreeFigureResult(
        organelle=figure_input.organelle,
        track_id=figure_input.track_id,
        sample_count=figure_input.sample_count,
        alignment_sites=figure_input.alignment_sites,
        model=figure_input.model,
        method=figure_input.method,
        treefile_path=figure_input.treefile_path,
        png_path=png_path,
        pdf_path=pdf_path,
        svg_path=svg_path,
        tip_count=tip_count,
        figure_width=width,
        figure_height=height,
    )


def write_tree_visualization_outputs(
    output_dir: Path,
    results: list[TreeFigureResult],
    run_label: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / labeled_output_name("tree_visualization_summary.tsv", run_label),
        [
            {
                "organelle": result.organelle,
                "track_id": result.track_id,
                "sample_count": str(result.sample_count),
                "tip_count": str(result.tip_count),
                "alignment_sites": str(result.alignment_sites),
                "method": result.method,
                "model": result.model,
                "treefile_path": result.treefile_path.as_posix(),
                "png_path": result.png_path.as_posix(),
                "pdf_path": result.pdf_path.as_posix(),
                "svg_path": result.svg_path.as_posix(),
                "figure_width": f"{result.figure_width:.2f}",
                "figure_height": f"{result.figure_height:.2f}",
            }
            for result in results
        ],
        [
            "organelle",
            "track_id",
            "sample_count",
            "tip_count",
            "alignment_sites",
            "method",
            "model",
            "treefile_path",
            "png_path",
            "pdf_path",
            "svg_path",
            "figure_width",
            "figure_height",
        ],
    )
    write_report(
        output_dir / labeled_output_name("tree_visualization_report.md", run_label),
        results=results,
        run_label=run_label,
    )


def write_report(path: Path, results: list[TreeFigureResult], run_label: str) -> None:
    label = run_label or "full"
    lines = [
        "# Step 14 Tree Visualizations",
        "",
        "This step renders the Step 11 IQ-TREE Newick trees into static figures",
        "for visual inspection, reporting, and downstream sharing.",
        "",
        "## Run",
        "",
        f"- Run label: `{label}`",
        "- Formats: PNG, PDF, SVG",
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
                f"- Samples in tree summary: {result.sample_count}",
                f"- Tips rendered: {result.tip_count}",
                f"- Alignment sites: {result.alignment_sites}",
                f"- Method: `{result.method}`",
                f"- Model: `{result.model}`",
                f"- PNG: `{result.png_path}`",
                f"- PDF: `{result.pdf_path}`",
                f"- SVG: `{result.svg_path}`",
                "",
            ]
        )
    path.write_text("\n".join(lines))


def run_tree_visualizations(
    tree_dir: Path = DEFAULT_TREE_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
) -> list[TreeFigureResult]:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/dudleya_matplotlib")
    inputs = read_tree_figure_inputs(tree_dir=tree_dir, run_label=run_label)
    results = [
        render_tree_figure(figure_input, output_dir=output_dir, run_label=run_label)
        for figure_input in inputs
    ]
    write_tree_visualization_outputs(output_dir, results, run_label=run_label)
    return results


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Step 14: render cpDNA/mtDNA phylogenetic tree figures."
    )
    parser.add_argument("--tree-dir", type=Path, default=DEFAULT_TREE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    results = run_tree_visualizations(
        tree_dir=args.tree_dir,
        output_dir=args.output_dir,
        run_label=args.run_label,
    )
    for result in results:
        print(
            f"{result.organelle}: rendered {result.tip_count} tips to "
            f"{result.png_path}, {result.pdf_path}, and {result.svg_path}"
        )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
