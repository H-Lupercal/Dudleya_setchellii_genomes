"""Build haploid cpDNA and mtDNA haplotype-network inputs and figures."""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path


BASES = frozenset("ACGT")
DEFAULT_SNP_ALIGNMENT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/10_snp_alignment"
)
DEFAULT_METADATA_PATH = Path(
    "dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/included_samples.tsv"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/21_haplotype_network"
)
DEFAULT_RENDERER_PATH = Path(
    "dudleya_organelle_alignment_pipeline/scripts/render_haplotype_network.R"
)
DEFAULT_RUN_LABEL = "primary"


class HaplotypeNetworkError(RuntimeError):
    """Raised when haplotype-network inputs or outputs are invalid."""


@dataclass(frozen=True)
class NetworkPaths:
    """All per-organelle Stage 21 paths sharing one output prefix."""

    prefix: Path
    input_fasta: Path
    site_table: Path
    metadata_table: Path
    popart_nexus: Path
    assignments: Path
    haplotype_summary: Path
    edges: Path
    layout: Path
    renderer_summary: Path
    png: Path
    pdf: Path
    svg: Path

    @property
    def renderer_outputs(self) -> tuple[Path, ...]:
        return (
            self.assignments,
            self.haplotype_summary,
            self.edges,
            self.layout,
            self.renderer_summary,
            self.png,
            self.pdf,
            self.svg,
        )


@dataclass(frozen=True)
class HaplotypeNetworkInput:
    organelle: str
    track_id: str
    sample_count: int
    source_site_count: int
    source_missing_bases: int
    alignment_fasta_path: Path
    site_table_path: Path


@dataclass(frozen=True)
class HaplotypeNetworkResult:
    organelle: str
    track_id: str
    sample_count: int
    source_site_count: int
    retained_site_count: int
    dropped_site_count: int
    haplotype_count: int
    edge_count: int
    species_group_count: int
    method: str
    paths: NetworkPaths


def network_paths(
    output_dir: Path,
    organelle: str,
    run_label: str,
) -> NetworkPaths:
    prefix = output_dir / f"{organelle}.{run_label}"
    return NetworkPaths(
        prefix=prefix,
        input_fasta=Path(f"{prefix}.haplotype_network_input.fa"),
        site_table=Path(f"{prefix}.haplotype_network_sites.tsv"),
        metadata_table=Path(f"{prefix}.haplotype_network_metadata.tsv"),
        popart_nexus=Path(f"{prefix}.popart.nex"),
        assignments=Path(f"{prefix}.haplotype_assignments.tsv"),
        haplotype_summary=Path(f"{prefix}.haplotype_summary.tsv"),
        edges=Path(f"{prefix}.haplotype_network_edges.tsv"),
        layout=Path(f"{prefix}.haplotype_network_layout.tsv"),
        renderer_summary=Path(
            f"{prefix}.haplotype_network_renderer_summary.tsv"
        ),
        png=Path(f"{prefix}.haplotype_network.png"),
        pdf=Path(f"{prefix}.haplotype_network.pdf"),
        svg=Path(f"{prefix}.haplotype_network.svg"),
    )


def _nexus_token(value: str) -> str:
    """Return a simple unquoted NEXUS token for project identifiers/traits."""

    return "_".join(value.split()).replace("'", "_")


def build_popart_nexus(
    records: list[tuple[str, str]],
    metadata: dict[str, dict[str, str]],
) -> str:
    """Build a PopART-readable DNA alignment with a species trait block."""

    if not records:
        raise HaplotypeNetworkError("No records supplied for PopART export")
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise HaplotypeNetworkError("PopART records have inconsistent lengths")
    character_count = lengths.pop()
    matrix = "\n".join(
        f"{_nexus_token(sample_id)} {sequence}"
        for sample_id, sequence in records
    )
    sample_groups = {
        sample_id: metadata[sample_id].get("species", "").strip()
        or "unresolved"
        for sample_id, _ in records
    }
    groups = list(dict.fromkeys(sample_groups.values()))
    traits = "\n".join(
        (
            f"{_nexus_token(sample_id)} "
            + ",".join(
                "1" if sample_groups[sample_id] == group else "0"
                for group in groups
            )
        )
        for sample_id, _ in records
    )
    return (
        "#NEXUS\n\n"
        "BEGIN DATA;\n"
        f"  DIMENSIONS NTAX={len(records)} NCHAR={character_count};\n"
        "  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n"
        "  MATRIX\n"
        f"{matrix}\n"
        "  ;\n"
        "END;\n\n"
        "BEGIN TRAITS;\n"
        f"  DIMENSIONS NTRAITS={len(groups)};\n"
        "  FORMAT LABELS=YES MISSING=? SEPARATOR=COMMA;\n"
        "  TRAITLABELS "
        + " ".join(_nexus_token(group) for group in groups)
        + ";\n"
        "  MATRIX\n"
        f"{traits}\n"
        "  ;\n"
        "END;\n"
    )


def build_renderer_command(
    rscript: Path,
    renderer: Path,
    fasta: Path,
    metadata: Path,
    prefix: Path,
    organelle: str,
) -> list[str]:
    return [
        str(rscript),
        str(renderer),
        str(fasta),
        str(metadata),
        str(prefix),
        organelle,
    ]


def validate_renderer_outputs(paths: NetworkPaths) -> None:
    """Reject missing or empty renderer products before reporting success."""

    for path in paths.renderer_outputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise HaplotypeNetworkError(
                f"Required renderer output is missing or empty: {path}"
            )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(
    path: Path,
    rows: list[dict[str, str]],
    fieldnames: list[str],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


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
        raise HaplotypeNetworkError(f"No FASTA records found in {path}")
    return records


def write_network_input_fasta(
    path: Path,
    records: list[tuple[str, str]],
) -> None:
    path.write_text(
        "".join(f">{sample_id}\n{sequence}\n" for sample_id, sequence in records)
    )


def write_network_site_table(
    path: Path,
    source_rows: list[dict[str, str]],
    kept_indexes: list[int],
) -> None:
    if not source_rows:
        raise HaplotypeNetworkError("No SNP site rows supplied")
    keep_set = set(kept_indexes)
    rows = [
        {
            **row,
            "source_alignment_index_0based": str(index),
            "network_status": (
                "retained" if index in keep_set else "dropped_missing"
            ),
        }
        for index, row in enumerate(source_rows)
    ]
    write_tsv(
        path,
        rows,
        [
            *source_rows[0].keys(),
            "source_alignment_index_0based",
            "network_status",
        ],
    )


def write_network_metadata(
    path: Path,
    records: list[tuple[str, str]],
    metadata: dict[str, dict[str, str]],
) -> None:
    rows: list[dict[str, str]] = []
    for sample_id, _ in records:
        row = metadata[sample_id]
        species = row.get("species", "").strip()
        rows.append(
            {
                "sample_id": sample_id,
                "species_group": species or "unresolved",
                "species": species,
                "popcode": row.get("popcode", ""),
                "population_name": row.get("population_name", ""),
                "naming_profile": row.get("naming_profile", ""),
            }
        )
    write_tsv(
        path,
        rows,
        [
            "sample_id",
            "species_group",
            "species",
            "popcode",
            "population_name",
            "naming_profile",
        ],
    )


def filter_complete_case_sites(
    records: list[tuple[str, str]],
) -> tuple[list[tuple[str, str]], list[int], list[int]]:
    """Remove every alignment column containing a non-ACGT state."""

    if not records:
        raise HaplotypeNetworkError("No FASTA records supplied")
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise HaplotypeNetworkError("FASTA records have inconsistent lengths")
    site_count = lengths.pop()
    kept = [
        index
        for index in range(site_count)
        if all(sequence[index] in BASES for _, sequence in records)
    ]
    keep_set = set(kept)
    dropped = [index for index in range(site_count) if index not in keep_set]
    if not kept:
        raise HaplotypeNetworkError("No complete-case SNP sites remain")
    filtered = [
        (sample_id, "".join(sequence[index] for index in kept))
        for sample_id, sequence in records
    ]
    if len({sequence for _, sequence in filtered}) < 2:
        raise HaplotypeNetworkError("Fewer than two haplotypes remain")
    return filtered, kept, dropped


def validate_sample_metadata(
    sample_ids: list[str],
    metadata: dict[str, dict[str, str]],
) -> None:
    """Require exact sample identity and order across the FASTA and metadata."""

    if sample_ids != list(metadata):
        raise HaplotypeNetworkError("FASTA and metadata sample IDs or order differ")


def _labeled_name(stem: str, run_label: str) -> str:
    return f"{run_label}.{stem}" if run_label else stem


def read_haplotype_network_inputs(
    snp_alignment_dir: Path,
    run_label: str,
) -> list[HaplotypeNetworkInput]:
    summary_path = snp_alignment_dir / _labeled_name(
        "snp_alignment_summary.tsv",
        run_label,
    )
    rows = read_tsv(summary_path)
    inputs: list[HaplotypeNetworkInput] = []
    for row in rows:
        fasta_path = Path(row["alignment_fasta_path"])
        site_table_path = Path(row["site_table_path"])
        if not fasta_path.is_file():
            raise HaplotypeNetworkError(f"Missing SNP alignment FASTA: {fasta_path}")
        if not site_table_path.is_file():
            raise HaplotypeNetworkError(
                f"Missing SNP alignment site table: {site_table_path}"
            )
        inputs.append(
            HaplotypeNetworkInput(
                organelle=row["organelle"],
                track_id=row["track_id"],
                sample_count=int(row["sample_count"]),
                source_site_count=int(row["alignment_sites"]),
                source_missing_bases=int(row["missing_bases"]),
                alignment_fasta_path=fasta_path,
                site_table_path=site_table_path,
            )
        )
    if [item.organelle for item in inputs] != ["cpDNA", "mtDNA"]:
        raise HaplotypeNetworkError(
            f"Expected cpDNA then mtDNA rows in {summary_path}"
        )
    return inputs


def read_sample_metadata(path: Path) -> dict[str, dict[str, str]]:
    rows = read_tsv(path)
    metadata: dict[str, dict[str, str]] = {}
    for row in rows:
        sample_id = row.get("sample_id", "")
        if not sample_id:
            raise HaplotypeNetworkError(f"Blank sample_id in {path}")
        if sample_id in metadata:
            raise HaplotypeNetworkError(f"Duplicate sample_id in {path}: {sample_id}")
        metadata[sample_id] = row
    if not metadata:
        raise HaplotypeNetworkError(f"No sample metadata rows found in {path}")
    return metadata


def read_renderer_result(
    network_input: HaplotypeNetworkInput,
    paths: NetworkPaths,
    retained_site_count: int,
    dropped_site_count: int,
) -> HaplotypeNetworkResult:
    summary_rows = read_tsv(paths.renderer_summary)
    if len(summary_rows) != 1:
        raise HaplotypeNetworkError(
            f"Expected one renderer summary row in {paths.renderer_summary}"
        )
    summary = summary_rows[0]
    if summary["organelle"] != network_input.organelle:
        raise HaplotypeNetworkError(
            f"Renderer organelle mismatch in {paths.renderer_summary}"
        )
    sample_count = int(summary["sample_count"])
    haplotype_count = int(summary["haplotype_count"])
    edge_count = int(summary["edge_count"])
    species_group_count = int(summary["species_group_count"])
    if sample_count != network_input.sample_count:
        raise HaplotypeNetworkError(
            f"Renderer sample count mismatch in {paths.renderer_summary}"
        )

    assignments = read_tsv(paths.assignments)
    assignment_ids = [row["sample_id"] for row in assignments]
    if len(assignments) != sample_count or len(set(assignment_ids)) != sample_count:
        raise HaplotypeNetworkError(
            f"Renderer assignments are incomplete or duplicated: {paths.assignments}"
        )
    haplotypes = read_tsv(paths.haplotype_summary)
    haplotype_ids = {row["haplotype_id"] for row in haplotypes}
    if len(haplotypes) != haplotype_count or len(haplotype_ids) != haplotype_count:
        raise HaplotypeNetworkError(
            f"Renderer haplotype summary count mismatch: {paths.haplotype_summary}"
        )
    edges = read_tsv(paths.edges)
    if len(edges) != edge_count or any(
        row["from_haplotype"] not in haplotype_ids
        or row["to_haplotype"] not in haplotype_ids
        for row in edges
    ):
        raise HaplotypeNetworkError(
            f"Renderer edge table is inconsistent: {paths.edges}"
        )
    if retained_site_count + dropped_site_count != network_input.source_site_count:
        raise HaplotypeNetworkError(
            f"Site-filter counts do not match {network_input.site_table_path}"
        )
    return HaplotypeNetworkResult(
        organelle=network_input.organelle,
        track_id=network_input.track_id,
        sample_count=sample_count,
        source_site_count=network_input.source_site_count,
        retained_site_count=retained_site_count,
        dropped_site_count=dropped_site_count,
        haplotype_count=haplotype_count,
        edge_count=edge_count,
        species_group_count=species_group_count,
        method="pegas::haploNet",
        paths=paths,
    )


def _result_row(result: HaplotypeNetworkResult) -> dict[str, str]:
    return {
        "organelle": result.organelle,
        "track_id": result.track_id,
        "sample_count": str(result.sample_count),
        "source_site_count": str(result.source_site_count),
        "retained_site_count": str(result.retained_site_count),
        "dropped_site_count": str(result.dropped_site_count),
        "haplotype_count": str(result.haplotype_count),
        "edge_count": str(result.edge_count),
        "species_group_count": str(result.species_group_count),
        "method": result.method,
        "input_fasta_path": result.paths.input_fasta.as_posix(),
        "site_table_path": result.paths.site_table.as_posix(),
        "metadata_table_path": result.paths.metadata_table.as_posix(),
        "assignments_path": result.paths.assignments.as_posix(),
        "haplotype_summary_path": result.paths.haplotype_summary.as_posix(),
        "edges_path": result.paths.edges.as_posix(),
        "layout_path": result.paths.layout.as_posix(),
        "png_path": result.paths.png.as_posix(),
        "pdf_path": result.paths.pdf.as_posix(),
        "svg_path": result.paths.svg.as_posix(),
        "popart_nexus_path": result.paths.popart_nexus.as_posix(),
    }


def write_haplotype_network_outputs(
    output_dir: Path,
    results: list[HaplotypeNetworkResult],
    commands: list[list[str]],
    run_label: str,
) -> None:
    summary_rows = [_result_row(result) for result in results]
    write_tsv(
        output_dir / _labeled_name("haplotype_network_summary.tsv", run_label),
        summary_rows,
        list(summary_rows[0]),
    )
    write_tsv(
        output_dir / _labeled_name("haplotype_network_commands.tsv", run_label),
        [
            {"organelle": result.organelle, "command": shlex.join(command)}
            for result, command in zip(results, commands, strict=True)
        ],
        ["organelle", "command"],
    )
    report_lines = [
        f"# Stage 21 Haplotype Networks ({run_label or 'unlabeled'})",
        "",
        "This stage applies complete-case filtering to the Stage 10 haploid SNP alignments, then builds separate cpDNA and mtDNA networks with `pegas::haploNet`.",
        "",
        "Nodes are haplotypes, node area is sample frequency, colored sectors are species groups, and edge labels are mutational steps. These are descriptions of sequence relationships and haplotype sharing, not ancestry proportions.",
        "",
    ]
    for result in results:
        report_lines.extend(
            [
                f"## {result.organelle}",
                "",
                f"- Samples: {result.sample_count}",
                (
                    f"- Complete SNP sites: {result.retained_site_count} of "
                    f"{result.source_site_count} (dropped {result.dropped_site_count})"
                ),
                f"- Haplotypes: {result.haplotype_count}",
                f"- Network edges: {result.edge_count}",
                f"- Figure: `{result.paths.png}`",
                f"- Assignments: `{result.paths.assignments}`",
                f"- PopART NEXUS: `{result.paths.popart_nexus}`",
                "",
            ]
        )
    report_lines.extend(
        [
            "## Interpretation limits",
            "",
            "Organelle sites are linked and represent a single nonrecombining lineage per organelle. The networks do not estimate nuclear population structure, admixture fractions, direction of gene flow, or the timing of shared ancestry.",
            "",
        ]
    )
    report_path = output_dir / _labeled_name(
        "haplotype_network_report.md",
        run_label,
    )
    report_path.write_text("\n".join(report_lines))


def run_haplotype_network_analysis(
    snp_alignment_dir: Path = DEFAULT_SNP_ALIGNMENT_DIR,
    metadata_path: Path = DEFAULT_METADATA_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
    rscript: Path = Path("Rscript"),
    renderer_path: Path = DEFAULT_RENDERER_PATH,
    runner=subprocess.run,
) -> list[HaplotypeNetworkResult]:
    inputs = read_haplotype_network_inputs(snp_alignment_dir, run_label)
    metadata = read_sample_metadata(metadata_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    results: list[HaplotypeNetworkResult] = []
    commands: list[list[str]] = []
    for network_input in inputs:
        records = read_fasta(network_input.alignment_fasta_path)
        validate_sample_metadata(
            [sample_id for sample_id, _ in records],
            metadata,
        )
        if len(records) != network_input.sample_count:
            raise HaplotypeNetworkError(
                f"FASTA sample count mismatch: {network_input.alignment_fasta_path}"
            )
        source_site_rows = read_tsv(network_input.site_table_path)
        if len(source_site_rows) != network_input.source_site_count:
            raise HaplotypeNetworkError(
                f"Site-table count mismatch: {network_input.site_table_path}"
            )
        filtered, kept, dropped = filter_complete_case_sites(records)
        paths = network_paths(output_dir, network_input.organelle, run_label)
        write_network_input_fasta(paths.input_fasta, filtered)
        write_network_site_table(paths.site_table, source_site_rows, kept)
        write_network_metadata(paths.metadata_table, filtered, metadata)
        paths.popart_nexus.write_text(build_popart_nexus(filtered, metadata))
        command = build_renderer_command(
            rscript,
            renderer_path,
            paths.input_fasta,
            paths.metadata_table,
            paths.prefix,
            network_input.organelle,
        )
        completed = runner(
            command,
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode:
            detail = completed.stderr or completed.stdout or "unknown R error"
            raise HaplotypeNetworkError(
                f"Haplotype renderer failed for {network_input.organelle}: {detail}"
            )
        validate_renderer_outputs(paths)
        results.append(
            read_renderer_result(
                network_input,
                paths,
                len(kept),
                len(dropped),
            )
        )
        commands.append(command)
    write_haplotype_network_outputs(output_dir, results, commands, run_label)
    return results


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build cpDNA and mtDNA haplotype networks with pegas."
    )
    parser.add_argument(
        "--snp-alignment-dir",
        type=Path,
        default=DEFAULT_SNP_ALIGNMENT_DIR,
    )
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    parser.add_argument("--rscript", type=Path, default=Path("Rscript"))
    parser.add_argument("--renderer", type=Path, default=DEFAULT_RENDERER_PATH)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    results = run_haplotype_network_analysis(
        snp_alignment_dir=args.snp_alignment_dir,
        metadata_path=args.metadata,
        output_dir=args.output_dir,
        run_label=args.run_label,
        rscript=args.rscript,
        renderer_path=args.renderer,
    )
    for result in results:
        print(
            f"{result.organelle}: {result.haplotype_count} haplotypes from "
            f"{result.sample_count} samples at {result.paths.png}"
        )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
