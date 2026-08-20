"""Concatenate matched cpDNA and mtDNA callable-consensus alignments."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

from dudleya_organelle_alignment_pipeline.variant_calling import labeled_output_name


DEFAULT_CPDNA_PATH = Path(
    "dudleya_organelle_alignment_pipeline/results/11_callable_consensus/"
    "cpDNA.primary.callable_consensus.fa"
)
DEFAULT_MTDNA_PATH = Path(
    "dudleya_organelle_alignment_pipeline/results/11_callable_consensus/"
    "mtDNA.primary.callable_consensus.fa"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/22_concatenated_consensus"
)
DEFAULT_RUN_LABEL = "primary"


class ConcatenatedConsensusError(RuntimeError):
    """Raised when cpDNA and mtDNA alignments cannot be safely concatenated."""


@dataclass(frozen=True)
class FastaAlignment:
    sample_names: tuple[str, ...]
    sequences: dict[str, str]
    sequence_length: int


@dataclass(frozen=True)
class ConcatenatedAlignment:
    sample_names: tuple[str, ...]
    sequences: dict[str, str]
    cpdna_length: int
    mtdna_length: int

    @property
    def combined_length(self) -> int:
        return self.cpdna_length + self.mtdna_length

    @property
    def mtdna_start(self) -> int:
        return self.cpdna_length + 1

    @property
    def missing_bases(self) -> int:
        return sum(sequence.count("N") for sequence in self.sequences.values())


@dataclass(frozen=True)
class ConcatenationResult:
    sample_count: int
    cpdna_length: int
    mtdna_length: int
    combined_length: int
    cpdna_missing_bases: int
    mtdna_missing_bases: int
    missing_bases: int
    cpdna_path: Path
    mtdna_path: Path
    fasta_path: Path

    @property
    def cpdna_end(self) -> int:
        return self.cpdna_length

    @property
    def mtdna_start(self) -> int:
        return self.cpdna_length + 1

    @property
    def mtdna_end(self) -> int:
        return self.combined_length


def read_fasta_alignment(path: Path) -> FastaAlignment:
    sample_names: list[str] = []
    sequence_parts: dict[str, list[str]] = {}
    sample_name: str | None = None
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                sample_name = line[1:].split()[0]
                if sample_name in sequence_parts:
                    raise ConcatenatedConsensusError(
                        f"Duplicate FASTA identifier {sample_name} in {path}"
                    )
                sample_names.append(sample_name)
                sequence_parts[sample_name] = []
            else:
                if sample_name is None:
                    raise ConcatenatedConsensusError(
                        f"FASTA sequence before header in {path}"
                    )
                sequence_parts[sample_name].append(line.upper())
    sequences = {
        name: "".join(sequence_parts[name])
        for name in sample_names
    }
    if not sequences:
        raise ConcatenatedConsensusError(f"No FASTA records found in {path}")
    empty_samples = [name for name, sequence in sequences.items() if not sequence]
    if empty_samples:
        raise ConcatenatedConsensusError(
            f"Empty FASTA sequence for {empty_samples[0]} in {path}"
        )
    sequence_lengths = {len(sequence) for sequence in sequences.values()}
    if len(sequence_lengths) != 1:
        raise ConcatenatedConsensusError(
            f"Inconsistent FASTA sequence lengths in {path}"
        )
    return FastaAlignment(
        sample_names=tuple(sample_names),
        sequences=sequences,
        sequence_length=sequence_lengths.pop(),
    )


def concatenate_consensus_alignments(
    cpdna_path: Path,
    mtdna_path: Path,
) -> ConcatenatedAlignment:
    cpdna = read_fasta_alignment(cpdna_path)
    mtdna = read_fasta_alignment(mtdna_path)
    cpdna_samples = set(cpdna.sample_names)
    mtdna_samples = set(mtdna.sample_names)
    if cpdna_samples != mtdna_samples:
        cpdna_only = sorted(cpdna_samples - mtdna_samples)
        mtdna_only = sorted(mtdna_samples - cpdna_samples)
        raise ConcatenatedConsensusError(
            "Sample identifier mismatch between cpDNA and mtDNA alignments: "
            f"cpDNA-only={cpdna_only}; mtDNA-only={mtdna_only}"
        )
    return ConcatenatedAlignment(
        sample_names=cpdna.sample_names,
        sequences={
            sample_name: cpdna.sequences[sample_name] + mtdna.sequences[sample_name]
            for sample_name in cpdna.sample_names
        },
        cpdna_length=cpdna.sequence_length,
        mtdna_length=mtdna.sequence_length,
    )


def write_fasta(path: Path, alignment: ConcatenatedAlignment) -> None:
    with path.open("w") as handle:
        for sample_name in alignment.sample_names:
            handle.write(f">{sample_name}\n")
            sequence = alignment.sequences[sample_name]
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(
    output_dir: Path,
    alignment: ConcatenatedAlignment,
    result: ConcatenationResult,
    run_label: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_fasta(result.fasta_path, alignment)
    compatibility_fields = [
        "organelle",
        "track_id",
        "sample_count",
        "consensus_length",
        "missing_bases",
        "alignment_fasta_path",
    ]
    write_tsv(
        output_dir / labeled_output_name("callable_consensus_summary.tsv", run_label),
        [
            {
                "organelle": "cpDNA_mtDNA",
                "track_id": "cpdna_then_mtdna",
                "sample_count": str(result.sample_count),
                "consensus_length": str(result.combined_length),
                "missing_bases": str(result.missing_bases),
                "alignment_fasta_path": result.fasta_path.as_posix(),
            }
        ],
        compatibility_fields,
    )
    detail_fields = [
        "organelle",
        "track_id",
        "sample_count",
        "cpDNA_length",
        "mtDNA_length",
        "combined_length",
        "cpDNA_start",
        "cpDNA_end",
        "mtDNA_start",
        "mtDNA_end",
        "cpDNA_missing_bases",
        "mtDNA_missing_bases",
        "combined_missing_bases",
        "cpDNA_fasta_path",
        "mtDNA_fasta_path",
        "alignment_fasta_path",
    ]
    write_tsv(
        output_dir
        / labeled_output_name("concatenated_consensus_summary.tsv", run_label),
        [
            {
                "organelle": "cpDNA_mtDNA",
                "track_id": "cpdna_then_mtdna",
                "sample_count": str(result.sample_count),
                "cpDNA_length": str(result.cpdna_length),
                "mtDNA_length": str(result.mtdna_length),
                "combined_length": str(result.combined_length),
                "cpDNA_start": "1",
                "cpDNA_end": str(result.cpdna_end),
                "mtDNA_start": str(result.mtdna_start),
                "mtDNA_end": str(result.mtdna_end),
                "cpDNA_missing_bases": str(result.cpdna_missing_bases),
                "mtDNA_missing_bases": str(result.mtdna_missing_bases),
                "combined_missing_bases": str(result.missing_bases),
                "cpDNA_fasta_path": result.cpdna_path.as_posix(),
                "mtDNA_fasta_path": result.mtdna_path.as_posix(),
                "alignment_fasta_path": result.fasta_path.as_posix(),
            }
        ],
        detail_fields,
    )
    report_path = output_dir / labeled_output_name(
        "concatenated_consensus_report.md",
        run_label,
    )
    report_path.write_text(
        "\n".join(
            [
                "# Concatenated cpDNA + mtDNA Consensus Alignment",
                "",
                "Each sample's mtDNA callable-consensus sequence was appended",
                "unchanged to the end of the same sample's cpDNA sequence.",
                "",
                f"- Samples: {result.sample_count}",
                f"- cpDNA positions: 1-{result.cpdna_end}",
                f"- mtDNA positions: {result.mtdna_start}-{result.mtdna_end}",
                f"- Combined length: {result.combined_length}",
                f"- Combined missing bases: {result.missing_bases}",
                f"- FASTA: `{result.fasta_path}`",
                "",
            ]
        )
    )


def run_concatenation(
    cpdna_path: Path = DEFAULT_CPDNA_PATH,
    mtdna_path: Path = DEFAULT_MTDNA_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
) -> ConcatenationResult:
    cpdna = read_fasta_alignment(cpdna_path)
    mtdna = read_fasta_alignment(mtdna_path)
    alignment = concatenate_consensus_alignments(cpdna_path, mtdna_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    label = f".{run_label}" if run_label else ""
    fasta_path = output_dir / f"cpDNA_mtDNA{label}.concatenated_consensus.fa"
    result = ConcatenationResult(
        sample_count=len(alignment.sample_names),
        cpdna_length=alignment.cpdna_length,
        mtdna_length=alignment.mtdna_length,
        combined_length=alignment.combined_length,
        cpdna_missing_bases=sum(
            sequence.count("N") for sequence in cpdna.sequences.values()
        ),
        mtdna_missing_bases=sum(
            sequence.count("N") for sequence in mtdna.sequences.values()
        ),
        missing_bases=alignment.missing_bases,
        cpdna_path=cpdna_path,
        mtdna_path=mtdna_path,
        fasta_path=fasta_path,
    )
    write_outputs(output_dir, alignment, result, run_label)
    return result


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Append each sample's mtDNA consensus to its cpDNA consensus."
    )
    parser.add_argument("--cpdna-path", type=Path, default=DEFAULT_CPDNA_PATH)
    parser.add_argument("--mtdna-path", type=Path, default=DEFAULT_MTDNA_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    result = run_concatenation(
        cpdna_path=args.cpdna_path,
        mtdna_path=args.mtdna_path,
        output_dir=args.output_dir,
        run_label=args.run_label,
    )
    print(f"Samples concatenated: {result.sample_count}")
    print(f"cpDNA length: {result.cpdna_length}")
    print(f"mtDNA length: {result.mtdna_length}")
    print(f"Combined length: {result.combined_length}")
    print(f"Output FASTA: {result.fasta_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
