"""Build haploid cpDNA and mtDNA haplotype-network inputs and figures."""

from __future__ import annotations

import csv
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
    traits = "\n".join(
        (
            f"{_nexus_token(sample_id)} "
            f"{_nexus_token(metadata[sample_id].get('species', '').strip() or 'unresolved')}"
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
        "  DIMENSIONS NTRAITS=1;\n"
        "  FORMAT LABELS=YES MISSING=? SEPARATOR=COMMA;\n"
        "  TRAITLABELS species_group;\n"
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
