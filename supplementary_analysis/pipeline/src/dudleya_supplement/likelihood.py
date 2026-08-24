"""IQ-TREE likelihood mapping and conditional NeighborNet execution."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

from .io import write_tsv
from .phylogeny import likelihood_decision, parse_split_nexus, supported_incompatible_pair


def build_likelihood_command(
    *,
    alignment: Path,
    model: str,
    quartets: int,
    seed: int,
    prefix: Path,
    threads: int,
) -> list[str]:
    """Build an IQ-TREE likelihood-mapping command with an explicit DNA type.

    Callable organelle alignments can contain enough missing data that IQ-TREE
    cannot reliably infer the alphabet even when the observed symbols are only
    A, C, G, T, and N.
    """
    return [
        "iqtree3",
        "-s",
        str(alignment),
        "-st",
        "DNA",
        "-m",
        model,
        "-lmap",
        str(quartets),
        "-n",
        "0",
        "-seed",
        str(seed),
        "-pre",
        str(prefix),
        "-nt",
        str(threads),
    ]


def parse_likelihood_report(path: Path) -> dict[str, float | int]:
    text = path.read_text()
    section = text.split("Quartet support of areas 1-7", 1)
    if len(section) != 2:
        raise ValueError(f"Missing seven-region likelihood mapping section: {path}")
    before_resolution = section[1].split("Quartet resolution per sequence", 1)[0]
    totals = [line for line in before_resolution.splitlines() if re.match(r"^\s+\d+\s+\d+\s+\(\s*\d", line)]
    if not totals:
        raise ValueError(f"Missing total seven-region row: {path}")
    values = re.findall(r"(\d+)\s+\(\s*([0-9.]+)\)", totals[-1])
    if len(values) != 7:
        raise ValueError(f"Expected seven likelihood regions, found {len(values)}: {path}")
    overall = re.search(
        r"Number of fully resolved\s+quartets.*?\(=([0-9.]+)%\).*?"
        r"Number of partly resolved quartets.*?\(=([0-9.]+)%\).*?"
        r"Number of unresolved\s+quartets.*?\(=([0-9.]+)%\)",
        text,
        re.DOTALL,
    )
    if overall is None:
        raise ValueError(f"Missing overall likelihood mapping fractions: {path}")
    result: dict[str, float | int] = {}
    for index, (count, percentage) in enumerate(values, 1):
        result[f"region_{index}_count"] = int(count)
        result[f"region_{index}_fraction"] = float(percentage) / 100.0
    result["resolved_fraction"] = float(overall.group(1)) / 100.0
    result["side_fraction"] = float(overall.group(2)) / 100.0
    result["center_fraction"] = float(overall.group(3)) / 100.0
    return result


def _run_outline(root: Path, run_id: str, organelle: str) -> list[Path]:
    outline = shutil.which("outline")
    if outline is None:
        raise RuntimeError("NeighborNet trigger met, but pinned SplitsPy 0.0.10 executable 'outline' is unavailable")
    source = root / f"canonical_publication/results/trees/publication-20260817/{organelle}.primary.mldist"
    output_dir = root / f"supplementary_analysis/results/phylogeny/{run_id}/neighbornet/{organelle}"
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix = output_dir / "distance_matrix.phy"
    matrix.write_bytes(source.read_bytes())
    image = output_dir / "neighbornet.png"
    nexus = output_dir / "neighbornet.splits.nex"
    graph = output_dir / "neighbornet.tgf"
    subprocess.run(
        [
            outline,
            "-o",
            str(image.relative_to(root)),
            "-n",
            str(nexus.relative_to(root)),
            "-t",
            str(graph.relative_to(root)),
            str(matrix.relative_to(root)),
        ],
        cwd=root,
        check=True,
    )
    version = output_dir / "SplitsPy.version.txt"
    version.write_text("SplitsPy==0.0.10\n")
    return [matrix, image, nexus, graph, version]


def run_likelihood_mapping(root: Path, run_id: str, config: dict[str, object]) -> list[Path]:
    settings = config["likelihood_mapping"]  # type: ignore[index]
    seeds = config["seeds"]  # type: ignore[index]
    output_dir = root / f"supplementary_analysis/results/phylogeny/{run_id}/likelihood_mapping"
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    rows: list[dict[str, object]] = []
    for organelle, model, seed in (
        ("chloroplast", settings["chloroplast_model"], seeds["cp_tree"]),  # type: ignore[index]
        ("mitochondria", settings["mitochondria_model"], seeds["mt_tree"]),  # type: ignore[index]
    ):
        alignment = root / f"canonical_publication/results/alignments/publication-20260817/{organelle}.callable_alignment.fa"
        work_prefix = root / f"supplementary_analysis/work/{run_id}/likelihood_mapping/{organelle}"
        work_prefix.parent.mkdir(parents=True, exist_ok=True)
        work_report = Path(f"{work_prefix}.iqtree")
        if not work_report.is_file():
            command = build_likelihood_command(
                alignment=alignment.relative_to(root),
                model=str(model),
                quartets=int(settings["quartets"]),  # type: ignore[index]
                seed=int(seed),
                prefix=work_prefix.relative_to(root),
                threads=8,
            )
            subprocess.run(command, cwd=root, check=True)
        copied = []
        for suffix in (".iqtree", ".lmap.svg", ".lmap.eps"):
            source = Path(f"{work_prefix}{suffix}")
            if source.is_file():
                destination = output_dir / f"{organelle}{suffix}"
                shutil.copyfile(source, destination)
                copied.append(destination)
        report = output_dir / f"{organelle}.iqtree"
        stats = parse_likelihood_report(report)
        taxa, splits = parse_split_nexus(root / f"canonical_publication/results/trees/publication-20260817/{organelle}.primary.splits.nex")
        conflict = supported_incompatible_pair(taxa, splits, minimum_frequency=float(settings["split_trigger"]))  # type: ignore[index]
        decision = likelihood_decision(
            center_fraction=float(stats["center_fraction"]),
            side_fraction=float(stats["side_fraction"]),
            has_supported_conflict=conflict is not None,
            center_limit=float(settings["center_limit"]),  # type: ignore[index]
            side_trigger=float(settings["side_trigger"]),  # type: ignore[index]
        )
        rows.append(
            {
                "organelle": organelle,
                "model": model,
                "seed": seed,
                "quartets": settings["quartets"],  # type: ignore[index]
                **stats,
                "bootstrap_split_taxa": len(taxa),
                "supported_incompatible_pair": "yes" if conflict else "no",
                "decision": decision,
            }
        )
        outputs.extend(copied)
        if decision == "RUN_NEIGHBORNET":
            outputs.extend(_run_outline(root, run_id, organelle))
    table = output_dir / "likelihood_mapping_summary.tsv"
    write_tsv(table, rows, list(rows[0]), root)
    outputs.append(table)
    return sorted(set(outputs))
