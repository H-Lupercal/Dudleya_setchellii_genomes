"""IQ-TREE likelihood mapping and conditional NeighborNet execution."""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

from Bio import SeqIO

from .io import write_tsv
from .phylogeny import likelihood_decision, parse_split_nexus, supported_incompatible_pair


def mask_restricted_sequences(records: dict[str, str], intervals: list[tuple[int, int]], *, expected_length: int) -> dict[str, str]:
    """Slice equal-length sequences using BED-style zero-based half-open intervals."""
    if not records:
        raise ValueError("Mask-restricted alignment requires at least one sequence")
    source_lengths = {len(sequence) for sequence in records.values()}
    if len(source_lengths) != 1:
        raise ValueError("Mask-restricted alignment sequences must have equal length")
    source_length = source_lengths.pop()
    if any(start < 0 or end <= start or end > source_length for start, end in intervals):
        raise ValueError("Mask interval lies outside the source alignment")
    observed_length = sum(end - start for start, end in intervals)
    if observed_length != expected_length:
        raise ValueError(f"Mask-restricted alignment length {observed_length}; expected {expected_length}")
    return {name: "".join(sequence[start:end] for start, end in intervals) for name, sequence in records.items()}


def run_command_logged(command: list[str], *, cwd: Path, log: Path) -> None:
    """Run a verbose command while preserving its complete screen output."""
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as handle:
        subprocess.run(command, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, check=True)


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
    values = re.findall(r"(\d+)\s+\(\s*([0-9.]+)\s*\)", totals[-1])
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


def parse_likelihood_diagnostics(path: Path) -> dict[str, int]:
    """Extract model-assumption warnings from IQ-TREE's run log."""
    text = path.read_text()
    alignment = re.search(r"Alignment has (\d+) sequences with", text)
    composition = re.search(r"(\d+) sequences failed composition chi2 test", text)
    ambiguity = re.search(r"WARNING: (\d+) sequences contain more than 50% gaps/ambiguity", text)
    if alignment is None or composition is None:
        raise ValueError(f"Missing alignment/composition diagnostics: {path}")
    return {
        "alignment_sequence_count": int(alignment.group(1)),
        "composition_failed_count": int(composition.group(1)),
        "over_50pct_ambiguity_count": int(ambiguity.group(1)) if ambiguity else 0,
    }


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
            run_command_logged(command, cwd=root, log=Path(f"{work_prefix}.screen.log"))
        copied = []
        for suffix in (".iqtree", ".log", ".lmap.svg", ".lmap.eps"):
            source = Path(f"{work_prefix}{suffix}")
            if source.is_file():
                destination = output_dir / f"{organelle}{suffix}"
                shutil.copyfile(source, destination)
                copied.append(destination)
        report = output_dir / f"{organelle}.iqtree"
        stats = parse_likelihood_report(report)
        diagnostics = parse_likelihood_diagnostics(Path(f"{work_prefix}.log"))
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
                **diagnostics,
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

    mask_path = root / "canonical_publication/references/masks/publication-20260817/mitochondria_high_confidence_sites.bed"
    intervals: list[tuple[int, int]] = []
    for line in mask_path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        record, start, end, *_ = line.split("\t")
        if record != "mitochondria":
            raise ValueError(f"Unexpected record in mitochondrial analysis mask: {record}")
        intervals.append((int(start), int(end)))
    expected_length = int(settings["mitochondria_mask_length"])  # type: ignore[index]
    source_alignment = root / "canonical_publication/results/alignments/publication-20260817/mitochondria.callable_alignment.fa"
    source_records = {record.id: str(record.seq).upper() for record in SeqIO.parse(source_alignment, "fasta")}
    restricted = mask_restricted_sequences(source_records, intervals, expected_length=expected_length)
    if len(restricted) != 271:
        raise RuntimeError(f"Expected 271 mitochondrial mask-restricted sequences, found {len(restricted)}")
    sensitivity_dir = root / f"supplementary_analysis/results/phylogeny/{run_id}/likelihood_mapping_sensitivity"
    sensitivity_dir.mkdir(parents=True, exist_ok=True)
    restricted_alignment = sensitivity_dir / "mitochondria.mask_restricted.fa"
    with restricted_alignment.open("w") as handle:
        for name, sequence in restricted.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    work_prefix = root / f"supplementary_analysis/work/{run_id}/likelihood_mapping/mitochondria_mask_restricted"
    work_prefix.parent.mkdir(parents=True, exist_ok=True)
    if not Path(f"{work_prefix}.iqtree").is_file():
        command = build_likelihood_command(
            alignment=restricted_alignment.relative_to(root),
            model=str(settings["mitochondria_model"]),  # type: ignore[index]
            quartets=int(settings["quartets"]),  # type: ignore[index]
            seed=int(seeds["mt_tree"]),  # type: ignore[index]
            prefix=work_prefix.relative_to(root),
            threads=8,
        )
        run_command_logged(command, cwd=root, log=Path(f"{work_prefix}.screen.log"))
    copied_sensitivity: list[Path] = []
    for suffix in (".iqtree", ".log", ".lmap.svg", ".lmap.eps"):
        source = Path(f"{work_prefix}{suffix}")
        if source.is_file():
            destination = sensitivity_dir / f"mitochondria.mask_restricted{suffix}"
            shutil.copyfile(source, destination)
            copied_sensitivity.append(destination)
    sensitivity_stats = parse_likelihood_report(sensitivity_dir / "mitochondria.mask_restricted.iqtree")
    sensitivity_diagnostics = parse_likelihood_diagnostics(Path(f"{work_prefix}.log"))
    primary_mt = next(row for row in rows if row["organelle"] == "mitochondria")
    resolved_change = float(sensitivity_stats["resolved_fraction"]) - float(primary_mt["resolved_fraction"])
    side_change = float(sensitivity_stats["side_fraction"]) - float(primary_mt["side_fraction"])
    center_change = float(sensitivity_stats["center_fraction"]) - float(primary_mt["center_fraction"])
    sensitivity_row: dict[str, object] = {
        "analysis": "mitochondria_mask_restricted",
        "organelle": "mitochondria",
        "model": settings["mitochondria_model"],  # type: ignore[index]
        "seed": seeds["mt_tree"],  # type: ignore[index]
        "quartets": settings["quartets"],  # type: ignore[index]
        "source_alignment_columns": len(next(iter(source_records.values()))),
        "restricted_alignment_columns": expected_length,
        **sensitivity_stats,
        **sensitivity_diagnostics,
        "resolved_fraction_change_from_primary": f"{resolved_change:.12g}",
        "side_fraction_change_from_primary": f"{side_change:.12g}",
        "center_fraction_change_from_primary": f"{center_change:.12g}",
        "network_role": "diagnostic_only_no_neighbornet",
    }
    sensitivity_table = sensitivity_dir / "likelihood_mapping_sensitivity.tsv"
    write_tsv(sensitivity_table, [sensitivity_row], list(sensitivity_row), root)
    outputs.extend([restricted_alignment, *copied_sensitivity, sensitivity_table])
    return sorted(set(outputs))
