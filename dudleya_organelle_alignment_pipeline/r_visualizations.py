"""Render additive R alternatives for existing analysis figures.

The module discovers tables and trees produced by earlier pipeline stages. It
does not recompute any biological analysis and never targets the existing
Matplotlib figure names.
"""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path


DEFAULT_PIPELINE_DIR = Path("dudleya_organelle_alignment_pipeline")
DEFAULT_RUN_LABEL = "primary"


class RVisualizationError(RuntimeError):
    """Raised when an additive R figure cannot be prepared safely."""


@dataclass(frozen=True)
class FigureJob:
    family: str
    organelle: str
    stage: str
    renderer_path: Path
    renderer_suffix: str
    source_paths: tuple[Path, ...]
    output_prefix: Path
    extra_args: tuple[str, ...] = ()


@dataclass(frozen=True)
class FigureResult:
    job: FigureJob
    outputs: tuple[Path, Path, Path]
    command: tuple[str, ...]
    status: str


def figure_outputs(
    output_prefix: Path,
    renderer_suffix: str,
) -> tuple[Path, Path, Path]:
    additive_prefix = Path(f"{output_prefix}.{renderer_suffix}")
    return tuple(
        Path(f"{additive_prefix}.{extension}")
        for extension in ("png", "pdf", "svg")
    )  # type: ignore[return-value]


def build_renderer_command(rscript: Path, job: FigureJob) -> list[str]:
    additive_prefix = Path(f"{job.output_prefix}.{job.renderer_suffix}")
    if job.family.startswith("admixture_"):
        mode = "structure" if job.family == "admixture_structure" else "cv"
        arguments = [mode, *(str(path) for path in job.source_paths), str(additive_prefix)]
        arguments.extend(job.extra_args)
    else:
        arguments = [*(str(path) for path in job.source_paths), str(additive_prefix)]
        arguments.extend(job.extra_args)
    return [str(rscript), str(job.renderer_path), *arguments]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _resolve_recorded_path(path_text: str, repo_root: Path) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else repo_root / path


def discover_figure_jobs(
    pipeline_dir: Path = DEFAULT_PIPELINE_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
) -> list[FigureJob]:
    pipeline_dir = pipeline_dir.resolve()
    repo_root = pipeline_dir.parent
    scripts_dir = pipeline_dir / "scripts"
    results_dir = pipeline_dir / "results"
    jobs: list[FigureJob] = []

    pca_dir = results_dir / "15_pca"
    for organelle in ("cpDNA", "mtDNA"):
        prefix = pca_dir / f"{organelle}.{run_label}.pca"
        jobs.append(
            FigureJob(
                family="pca",
                organelle=organelle,
                stage="15_pca",
                renderer_path=scripts_dir / "render_pca_ggplot.R",
                renderer_suffix="r_ggplot",
                source_paths=(
                    Path(f"{prefix}.coordinates.tsv"),
                    Path(f"{prefix}.variance.tsv"),
                ),
                output_prefix=prefix,
                extra_args=(organelle,),
            )
        )

    for stage in ("16_admixture", "18_admixture_replicates"):
        stage_dir = results_dir / stage
        summary_path = stage_dir / f"{run_label}.admixture_summary.tsv"
        for organelle in ("cpDNA", "mtDNA"):
            q_tables = sorted(stage_dir.glob(f"{organelle}.{run_label}.bestK*.q.tsv"))
            if len(q_tables) != 1:
                raise RVisualizationError(
                    f"Expected one {organelle} Q table in {stage_dir}, found {len(q_tables)}"
                )
            q_table = q_tables[0]
            structure_name = q_table.name.removesuffix(".q.tsv") + ".structure"
            jobs.append(
                FigureJob(
                    family="admixture_structure",
                    organelle=organelle,
                    stage=stage,
                    renderer_path=scripts_dir / "render_admixture_ggplot.R",
                    renderer_suffix="r_ggplot",
                    source_paths=(q_table,),
                    output_prefix=stage_dir / structure_name,
                    extra_args=(organelle,),
                )
            )
            jobs.append(
                FigureJob(
                    family="admixture_cv",
                    organelle=organelle,
                    stage=stage,
                    renderer_path=scripts_dir / "render_admixture_ggplot.R",
                    renderer_suffix="r_ggplot",
                    source_paths=(summary_path,),
                    output_prefix=stage_dir / f"{organelle}.{run_label}.admixture_cv",
                    extra_args=(organelle,),
                )
            )

    tree_stages = (
        ("14_tree_visualization", "initial", "0"),
        ("20_bootstrap_tree_visualization", "bootstrap", "1000"),
    )
    metadata_path = results_dir / "07_downstream_sample_set" / "included_samples.tsv"
    for stage, mode, replicates in tree_stages:
        stage_dir = results_dir / stage
        summary_path = stage_dir / f"{run_label}.tree_visualization_summary.tsv"
        rows = _read_tsv(summary_path)
        for row in rows:
            organelle = row["organelle"]
            jobs.append(
                FigureJob(
                    family="tree",
                    organelle=organelle,
                    stage=stage,
                    renderer_path=scripts_dir / "render_tree_ggtree.R",
                    renderer_suffix="r_ggtree",
                    source_paths=(
                        _resolve_recorded_path(row["treefile_path"], repo_root),
                        metadata_path,
                    ),
                    output_prefix=stage_dir / f"{organelle}.{run_label}.iqtree_ml_tree",
                    extra_args=(organelle, mode, replicates),
                )
            )

    return jobs


def _write_stage_records(
    results: list[FigureResult],
    run_label: str,
) -> None:
    by_stage: dict[str, list[FigureResult]] = {}
    for result in results:
        by_stage.setdefault(result.job.stage, []).append(result)

    for stage_results in by_stage.values():
        stage_dir = stage_results[0].job.output_prefix.parent
        command_path = stage_dir / f"{run_label}.r_visualization_commands.tsv"
        with command_path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "family",
                    "organelle",
                    "status",
                    "source_paths",
                    "png_path",
                    "pdf_path",
                    "svg_path",
                    "command",
                ],
                delimiter="\t",
            )
            writer.writeheader()
            for result in stage_results:
                writer.writerow(
                    {
                        "family": result.job.family,
                        "organelle": result.job.organelle,
                        "status": result.status,
                        "source_paths": ";".join(
                            path.as_posix() for path in result.job.source_paths
                        ),
                        "png_path": result.outputs[0].as_posix(),
                        "pdf_path": result.outputs[1].as_posix(),
                        "svg_path": result.outputs[2].as_posix(),
                        "command": shlex.join(result.command),
                    }
                )

        report_path = stage_dir / f"{run_label}.r_visualization_report.md"
        lines = [
            f"# Additive R Visualizations: {stage_results[0].job.stage}",
            "",
            "These figures visualize existing pipeline results with R. They do not",
            "rerun or replace the biological analyses or the original figures.",
            "",
            "## Legend and key rules",
            "",
            "- PCA and trees use the fixed species-group palette; unresolved metadata is gray.",
            "- ADMIXTURE colors are inferred clusters with arbitrary labels, not named populations.",
            "- CV plots mark the lowest mean error and state that lower is better.",
            "- Bootstrap-tree internal values are UFBoot support percentages from 1,000 replicates.",
            "",
            "## Added figures",
            "",
        ]
        for result in stage_results:
            lines.extend(
                [
                    f"### {result.job.organelle} {result.job.family}",
                    "",
                    f"- Status: `{result.status}`",
                    f"- PNG: `{result.outputs[0]}`",
                    f"- PDF: `{result.outputs[1]}`",
                    f"- SVG: `{result.outputs[2]}`",
                    "",
                ]
            )
        report_path.write_text("\n".join(lines))


def run_r_visualizations(
    pipeline_dir: Path = DEFAULT_PIPELINE_DIR,
    rscript: Path = Path("Rscript"),
    run_label: str = DEFAULT_RUN_LABEL,
    stages: list[str] | None = None,
    force: bool = False,
) -> list[FigureResult]:
    pipeline_dir = pipeline_dir.resolve()
    repo_root = pipeline_dir.parent
    jobs = discover_figure_jobs(pipeline_dir=pipeline_dir, run_label=run_label)
    if stages:
        selected = set(stages)
        jobs = [job for job in jobs if job.stage in selected]
    if not jobs:
        raise RVisualizationError("No R visualization jobs matched the selection")

    results: list[FigureResult] = []
    for job in jobs:
        missing = [path for path in (job.renderer_path, *job.source_paths) if not path.exists()]
        if missing:
            raise RVisualizationError(
                "Missing R visualization input(s): "
                + ", ".join(path.as_posix() for path in missing)
            )
        outputs = figure_outputs(job.output_prefix, job.renderer_suffix)
        existing = [path.exists() for path in outputs]
        command = build_renderer_command(rscript, job)
        if all(existing) and not force:
            status = "outputs already present; pass --force to regenerate"
        else:
            if any(existing) and not force:
                raise RVisualizationError(
                    f"Partial additive output set exists for {job.output_prefix}; "
                    "pass --force to regenerate all formats"
                )
            completed = subprocess.run(
                command,
                cwd=repo_root,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            if completed.returncode:
                raise RVisualizationError(
                    f"R renderer failed for {job.stage} {job.organelle}:\n"
                    f"{completed.stdout}"
                )
            invalid = [path for path in outputs if not path.exists() or path.stat().st_size == 0]
            if invalid:
                raise RVisualizationError(
                    "R renderer did not write non-empty output(s): "
                    + ", ".join(path.as_posix() for path in invalid)
                )
            status = "rendered"
        results.append(
            FigureResult(
                job=job,
                outputs=outputs,
                command=tuple(command),
                status=status,
            )
        )

    _write_stage_records(results, run_label=run_label)
    return results


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Render additive R alternatives for existing pipeline figures."
    )
    parser.add_argument("--pipeline-dir", type=Path, default=DEFAULT_PIPELINE_DIR)
    parser.add_argument("--rscript", type=Path, default=Path("Rscript"))
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    parser.add_argument(
        "--stages",
        nargs="+",
        choices=[
            "14_tree_visualization",
            "15_pca",
            "16_admixture",
            "18_admixture_replicates",
            "20_bootstrap_tree_visualization",
        ],
    )
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    results = run_r_visualizations(
        pipeline_dir=args.pipeline_dir,
        rscript=args.rscript,
        run_label=args.run_label,
        stages=args.stages,
        force=args.force,
    )
    for result in results:
        print(
            f"{result.job.stage} {result.job.organelle} {result.job.family}: "
            f"{result.status} -> {result.outputs[0]}"
        )
    print(f"Added or verified {len(results)} R figure sets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
