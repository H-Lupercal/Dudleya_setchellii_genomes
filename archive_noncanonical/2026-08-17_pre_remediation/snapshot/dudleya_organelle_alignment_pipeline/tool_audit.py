"""Audit local bioinformatics tools needed for the organelle popgen workflow."""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


DEFAULT_OUTPUT_DIR = Path("dudleya_organelle_alignment_pipeline/results/13_tool_audit")
DEFAULT_AUDIT_LABEL = "primary"


@dataclass(frozen=True)
class ToolSpec:
    tool_id: str
    executables: tuple[str, ...]
    necessity: str
    required_for: str
    version_args: tuple[str, ...]


@dataclass(frozen=True)
class ToolResult:
    tool_id: str
    executable: str
    necessity: str
    required_for: str
    status: str
    path: str
    version: str
    note: str


@dataclass(frozen=True)
class AuditSummary:
    total_tools: int
    found_tools: int
    missing_required_current: list[str]
    missing_required_remaining: list[str]
    missing_recommended: list[str]

    @property
    def ready_for_current_pipeline(self) -> bool:
        return not self.missing_required_current

    @property
    def ready_for_remaining_goal(self) -> bool:
        return not self.missing_required_current and not self.missing_required_remaining


Resolver = Callable[[str], str | None]
Runner = Callable[[list[str]], str]


TOOL_SPECS = (
    ToolSpec(
        "python3",
        ("python3",),
        "required_current",
        "pipeline scripts, tests, PCA/Fst helper code",
        ("--version",),
    ),
    ToolSpec(
        "bwa",
        ("bwa",),
        "required_current",
        "read mapping to combined cpDNA/mtDNA reference",
        (),
    ),
    ToolSpec(
        "samtools",
        ("samtools",),
        "required_current",
        "BAM sorting/indexing/depth/QC",
        ("--version",),
    ),
    ToolSpec(
        "bcftools",
        ("bcftools",),
        "required_current",
        "haploid variant calling, filtering, and VCF indexing",
        ("--version",),
    ),
    ToolSpec(
        "fastp",
        ("fastp",),
        "required_current",
        "read QC/trimming if rerunning raw-read QC",
        ("--version",),
    ),
    ToolSpec(
        "fastqc",
        ("fastqc",),
        "required_current",
        "read QC reports",
        ("--version",),
    ),
    ToolSpec(
        "multiqc",
        ("multiqc",),
        "required_current",
        "aggregate QC reports",
        ("--version",),
    ),
    ToolSpec(
        "iqtree",
        ("iqtree", "iqtree2"),
        "required_current",
        "maximum-likelihood cpDNA/mtDNA phylogenetic trees",
        ("--version",),
    ),
    ToolSpec(
        "FastTree",
        ("FastTree", "fasttree"),
        "recommended_remaining",
        "quick approximate ML tree checks",
        ("-help",),
    ),
    ToolSpec(
        "plink",
        ("plink", "plink2"),
        "required_remaining",
        "PCA matrix handling and ADMIXTURE input preparation",
        ("--version",),
    ),
    ToolSpec(
        "admixture",
        ("admixture",),
        "required_remaining",
        "structure/admixture-style clustering and empirical K selection",
        ("--version",),
    ),
    ToolSpec(
        "vcftools",
        ("vcftools",),
        "recommended_remaining",
        "VCF-based population summary checks",
        ("--version",),
    ),
    ToolSpec(
        "bedtools",
        ("bedtools",),
        "recommended_remaining",
        "interval QC and mask cross-checking",
        ("--version",),
    ),
    ToolSpec(
        "Rscript",
        ("Rscript",),
        "required_remaining",
        "PCA/tree/admixture/Fst plotting outputs",
        ("--version",),
    ),
    ToolSpec(
        "snakemake",
        ("snakemake",),
        "recommended_remaining",
        "optional integration with Snakemake orchestration",
        ("--version",),
    ),
    ToolSpec(
        "python_matplotlib",
        ("python3",),
        "required_remaining",
        "PCA scatterplots, tree rendering, and static PNG/PDF figures",
        ("-c", "import matplotlib; print(matplotlib.__version__)"),
    ),
    ToolSpec(
        "python_pandas",
        ("python3",),
        "required_remaining",
        "figure-ready metadata tables and plotting data frames",
        ("-c", "import pandas; print(pandas.__version__)"),
    ),
    ToolSpec(
        "python_sklearn",
        ("python3",),
        "required_remaining",
        "PCA calculation and variance summaries",
        ("-c", "import sklearn; print(sklearn.__version__)"),
    ),
    ToolSpec(
        "python_biopython",
        ("python3",),
        "required_remaining",
        "Newick tree parsing and scripted tree rendering",
        ("-c", "import Bio; print(Bio.__version__)"),
    ),
    ToolSpec(
        "python_seaborn",
        ("python3",),
        "recommended_remaining",
        "polished statistical plots",
        ("-c", "import seaborn; print(seaborn.__version__)"),
    ),
    ToolSpec(
        "python_ete3",
        ("python3",),
        "recommended_remaining",
        "alternate Newick tree visualization",
        ("-c", "import ete3; print(ete3.__version__)"),
    ),
    ToolSpec(
        "r_ggplot2",
        ("Rscript",),
        "required_remaining",
        "R-based PCA/admixture/Fst plots",
        ("-e", "cat(as.character(packageVersion('ggplot2')))"),
    ),
    ToolSpec(
        "r_ape",
        ("Rscript",),
        "required_remaining",
        "R-based phylogenetic tree parsing/rendering",
        ("-e", "cat(as.character(packageVersion('ape')))"),
    ),
    ToolSpec(
        "r_ggtree",
        ("Rscript",),
        "required_current",
        "additive R phylogenetic tree figures with support and species keys",
        ("-e", "cat(as.character(packageVersion('ggtree')))"),
    ),
    ToolSpec(
        "r_pegas",
        ("Rscript",),
        "required_current",
        "haploid cpDNA/mtDNA haplotype networks",
        ("-e", "cat(as.character(packageVersion('pegas')))"),
    ),
    ToolSpec(
        "r_patchwork",
        ("Rscript",),
        "recommended_remaining",
        "combining multiple R figure panels",
        ("-e", "cat(as.character(packageVersion('patchwork')))"),
    ),
)


def run_version_command(command: list[str]) -> str:
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", "/tmp/dudleya_matplotlib")
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        env=env,
    )
    if completed.returncode:
        output = first_version_line(completed.stdout.strip())
        raise RuntimeError(output or f"command exited {completed.returncode}")
    return completed.stdout.strip()


def first_version_line(output: str) -> str:
    for line in output.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return ""


def missing_note(necessity: str) -> str:
    if necessity == "required_current":
        return "Install before rerunning the completed mapping/QC/variant/tree pipeline."
    if necessity == "required_remaining":
        return "Install before continuing to the remaining planned analyses."
    return "Recommended, but not strictly blocking for the immediate next scripted step."


def check_tool(
    spec: ToolSpec,
    resolver: Resolver = shutil.which,
    runner: Runner = run_version_command,
) -> ToolResult:
    for executable in spec.executables:
        path = resolver(executable)
        if path:
            command = [path, *spec.version_args]
            try:
                version = first_version_line(runner(command)) if spec.version_args else executable
            except Exception as exc:
                return ToolResult(
                    tool_id=spec.tool_id,
                    executable=executable,
                    necessity=spec.necessity,
                    required_for=spec.required_for,
                    status="MISSING",
                    path=path,
                    version="",
                    note=f"Executable found, but version/import check failed: {exc}",
                )
            return ToolResult(
                tool_id=spec.tool_id,
                executable=executable,
                necessity=spec.necessity,
                required_for=spec.required_for,
                status="FOUND",
                path=path,
                version=version or "version_not_reported",
                note="Available on PATH.",
            )
    return ToolResult(
        tool_id=spec.tool_id,
        executable=";".join(spec.executables),
        necessity=spec.necessity,
        required_for=spec.required_for,
        status="MISSING",
        path="",
        version="",
        note=missing_note(spec.necessity),
    )


def run_tool_audit(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    audit_label: str = DEFAULT_AUDIT_LABEL,
) -> list[ToolResult]:
    results = [check_tool(spec) for spec in TOOL_SPECS]
    write_tool_audit_outputs(output_dir, results, audit_label)
    return results


def summarize_audit(results: list[ToolResult]) -> AuditSummary:
    return AuditSummary(
        total_tools=len(results),
        found_tools=sum(1 for result in results if result.status == "FOUND"),
        missing_required_current=[
            result.tool_id
            for result in results
            if result.status == "MISSING" and result.necessity == "required_current"
        ],
        missing_required_remaining=[
            result.tool_id
            for result in results
            if result.status == "MISSING" and result.necessity == "required_remaining"
        ],
        missing_recommended=[
            result.tool_id
            for result in results
            if result.status == "MISSING" and result.necessity == "recommended_remaining"
        ],
    )


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_tool_audit_outputs(
    output_dir: Path,
    results: list[ToolResult],
    audit_label: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / f"{audit_label}.tool_audit.tsv",
        [
            {
                "tool_id": result.tool_id,
                "executable": result.executable,
                "necessity": result.necessity,
                "required_for": result.required_for,
                "status": result.status,
                "path": result.path,
                "version": result.version,
                "note": result.note,
            }
            for result in results
        ],
        [
            "tool_id",
            "executable",
            "necessity",
            "required_for",
            "status",
            "path",
            "version",
            "note",
        ],
    )
    write_report(output_dir / f"{audit_label}.tool_audit_report.md", results, audit_label)


def write_report(path: Path, results: list[ToolResult], audit_label: str) -> None:
    summary = summarize_audit(results)
    lines = [
        "# Bioinformatics Tool Audit",
        "",
        "This report records whether the local tools needed for the Dudleya",
        "cpDNA/mtDNA workflow are installed and visible on `PATH`.",
        "",
        "## Summary",
        "",
        f"- Audit label: `{audit_label}`",
        f"- Tools checked: {summary.total_tools}",
        f"- Tools found: {summary.found_tools}",
        f"- Missing required current-pipeline tools: {', '.join(summary.missing_required_current) or 'none'}",
        f"- Missing required remaining-goal tools: {', '.join(summary.missing_required_remaining) or 'none'}",
        f"- Missing recommended tools: {', '.join(summary.missing_recommended) or 'none'}",
        "",
        "## Tool Checks",
        "",
        "| Tool | Necessity | Status | Version | Path | Use |",
        "|---|---|---|---|---|---|",
    ]
    for result in results:
        lines.append(
            "| "
            + " | ".join(
                [
                    result.tool_id,
                    result.necessity,
                    result.status,
                    result.version or "",
                    result.path or "",
                    result.required_for,
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            (
                "- The completed mapping, variant, consensus, and first-pass tree steps "
                "can be reproduced with the required current-pipeline tools."
                if summary.ready_for_current_pipeline
                else "- The completed pipeline is not fully reproducible until the missing current-pipeline tools are installed."
            ),
            (
                "- The remaining planned analyses have their required external tools available."
                if summary.ready_for_remaining_goal
                else "- Do not continue to the remaining planned analyses until missing required remaining-goal tools are installed or replaced by documented in-repo implementations."
            ),
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Audit Dudleya organelle workflow tools.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--audit-label", default=DEFAULT_AUDIT_LABEL)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    results = run_tool_audit(output_dir=args.output_dir, audit_label=args.audit_label)
    summary = summarize_audit(results)
    print(f"Tools checked: {summary.total_tools}")
    print(f"Tools found: {summary.found_tools}")
    print(
        "Missing required remaining-goal tools: "
        + (", ".join(summary.missing_required_remaining) or "none")
    )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
