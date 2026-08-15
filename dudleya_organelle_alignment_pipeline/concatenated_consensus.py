"""Concatenate matched cpDNA and mtDNA callable-consensus alignments."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


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
