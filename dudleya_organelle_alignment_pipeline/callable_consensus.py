"""Build full callable-site cpDNA/mtDNA consensus FASTA alignments.

This is Step 10 of the pipeline. It consumes Step 7 raw haploid variants,
Step 8 filtered SNPs, Step 4 population-genetic BED tracks, and Step 5 depth
files to produce one all-sample consensus FASTA per organelle.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from dudleya_organelle_alignment_pipeline.pilot_alignment import safe_sample_name
from dudleya_organelle_alignment_pipeline.variant_calling import labeled_output_name


DEFAULT_REFERENCE = Path(
    "dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa"
)
DEFAULT_SAMPLE_TABLE = Path(
    "dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/"
    "included_samples.tsv"
)
DEFAULT_TRACK_TABLE = Path(
    "dudleya_organelle_alignment_pipeline/results/05_analysis_masks/analysis_tracks.tsv"
)
DEFAULT_DEPTH_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/qc"
)
DEFAULT_VARIANT_CALLING_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/08_variant_calling"
)
DEFAULT_VARIANT_FILTERING_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/09_variant_filtering"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/11_callable_consensus"
)
DEFAULT_RUN_LABEL = "primary"
DEFAULT_MIN_DEPTH = 1


class CallableConsensusError(RuntimeError):
    """Raised when Step 10 cannot safely build consensus alignments."""


@dataclass(frozen=True)
class BedInterval:
    chrom: str
    start_0based: int
    end_0based: int
    name: str

    @property
    def positions_1based(self) -> range:
        return range(self.start_0based + 1, self.end_0based + 1)

    @property
    def length(self) -> int:
        return self.end_0based - self.start_0based


@dataclass(frozen=True)
class VcfVariant:
    chrom: str
    position: int
    ref: str
    alt: str
    genotypes: dict[str, str]


@dataclass(frozen=True)
class CallableAlignment:
    sample_names: list[str]
    sequences: dict[str, str]
    sites: list[dict[str, str]]
    filtered_variant_sites: int
    masked_failed_variant_sites: int

    @property
    def consensus_length(self) -> int:
        return len(self.sites)

    @property
    def missing_bases(self) -> int:
        return sum(sequence.count("N") for sequence in self.sequences.values())


@dataclass(frozen=True)
class ConsensusInput:
    organelle: str
    track_id: str
    bed_path: Path
    raw_records: int
    filtered_records: int
    raw_vcf_path: Path
    filtered_vcf_path: Path

    def to_result(
        self,
        sample_count: int,
        consensus_length: int,
        filtered_variant_sites: int,
        masked_failed_variant_sites: int,
        missing_bases: int,
        fasta_path: Path,
        site_table_path: Path,
    ) -> "ConsensusResult":
        return ConsensusResult(
            organelle=self.organelle,
            track_id=self.track_id,
            sample_count=sample_count,
            consensus_length=consensus_length,
            raw_records=self.raw_records,
            filtered_records=self.filtered_records,
            filtered_variant_sites=filtered_variant_sites,
            masked_failed_variant_sites=masked_failed_variant_sites,
            missing_bases=missing_bases,
            bed_path=self.bed_path,
            raw_vcf_path=self.raw_vcf_path,
            filtered_vcf_path=self.filtered_vcf_path,
            fasta_path=fasta_path,
            site_table_path=site_table_path,
        )


@dataclass(frozen=True)
class ConsensusResult:
    organelle: str
    track_id: str
    sample_count: int
    consensus_length: int
    raw_records: int
    filtered_records: int
    filtered_variant_sites: int
    masked_failed_variant_sites: int
    missing_bases: int
    bed_path: Path
    raw_vcf_path: Path
    filtered_vcf_path: Path
    fasta_path: Path
    site_table_path: Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
                continue
            if name is None:
                raise CallableConsensusError(f"FASTA sequence before header in {path}")
            records[name].append(line.upper())
    if not records:
        raise CallableConsensusError(f"No FASTA records found in {path}")
    return {record: "".join(parts) for record, parts in records.items()}


def read_bed(path: Path) -> list[BedInterval]:
    intervals: list[BedInterval] = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if not fields or fields[0].startswith("#"):
                continue
            if len(fields) < 3:
                raise CallableConsensusError(f"Invalid BED row in {path}: {fields}")
            start = int(fields[1])
            end = int(fields[2])
            if start < 0 or end <= start:
                raise CallableConsensusError(f"Invalid BED interval in {path}: {fields}")
            intervals.append(
                BedInterval(
                    chrom=fields[0],
                    start_0based=start,
                    end_0based=end,
                    name=fields[3] if len(fields) > 3 else f"region_{len(intervals) + 1}",
                )
            )
    if not intervals:
        raise CallableConsensusError(f"No BED intervals found in {path}")
    return intervals


def read_sample_names(sample_table: Path, organelle: str) -> list[str]:
    use_column = "downstream_cpDNA_use" if organelle == "cpDNA" else "downstream_mtDNA_use"
    sample_names = [
        row["sample_id"]
        for row in read_tsv(sample_table)
        if row.get(use_column) == "include"
    ]
    if not sample_names:
        raise CallableConsensusError(
            f"No included {organelle} samples found in {sample_table}"
        )
    return sample_names


def read_vcf_variants(path: Path) -> tuple[list[str], list[VcfVariant]]:
    sample_names: list[str] | None = None
    variants: list[VcfVariant] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_names = line.split("\t")[9:]
                continue
            if sample_names is None:
                raise CallableConsensusError(f"VCF header missing #CHROM line: {path}")
            fields = line.split("\t")
            chrom, position, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            variants.append(
                VcfVariant(
                    chrom=chrom,
                    position=position,
                    ref=ref.upper(),
                    alt=alt.upper(),
                    genotypes=dict(zip(sample_names, fields[9:], strict=True)),
                )
            )
    if sample_names is None:
        raise CallableConsensusError(f"VCF header missing #CHROM line: {path}")
    return sample_names, variants


def genotype_to_base(genotype_field: str, ref: str, alt: str) -> str:
    genotype = genotype_field.split(":", 1)[0]
    allele = genotype.replace("|", "/").split("/", 1)[0]
    if allele == "0":
        return ref.upper()
    if allele == "1":
        return alt.upper()
    return "N"


def build_track_template(
    reference: dict[str, str],
    intervals: list[BedInterval],
) -> tuple[list[str], dict[tuple[str, int], int], list[dict[str, str]]]:
    bases: list[str] = []
    coordinate_to_index: dict[tuple[str, int], int] = {}
    sites: list[dict[str, str]] = []
    for interval in intervals:
        if interval.chrom not in reference:
            raise CallableConsensusError(
                f"BED record {interval.chrom} is absent from the reference"
            )
        record = reference[interval.chrom]
        if interval.end_0based > len(record):
            raise CallableConsensusError(
                f"BED interval {interval.chrom}:{interval.start_0based}-{interval.end_0based} "
                f"exceeds reference length {len(record)}"
            )
        for position in interval.positions_1based:
            base = record[position - 1].upper()
            coordinate_to_index[(interval.chrom, position)] = len(bases)
            bases.append(base)
            sites.append(
                {
                    "site_index": str(len(sites) + 1),
                    "chrom": interval.chrom,
                    "position": str(position),
                    "reference_base": base,
                }
            )
    return bases, coordinate_to_index, sites


def read_covered_indexes(
    depth_path: Path,
    coordinate_to_index: dict[tuple[str, int], int],
    min_depth: int,
) -> set[int]:
    covered: set[int] = set()
    with depth_path.open() as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            index = coordinate_to_index.get((fields[0], int(fields[1])))
            if index is not None and int(fields[2]) >= min_depth:
                covered.add(index)
    return covered


def depth_path_for_sample(depth_dir: Path, sample_name: str) -> Path:
    return depth_dir / f"{safe_sample_name(sample_name)}.depth.tsv"


def build_callable_consensus(
    reference_path: Path,
    bed_path: Path,
    sample_table: Path,
    depth_dir: Path,
    raw_vcf_path: Path,
    filtered_vcf_path: Path,
    min_depth: int,
    organelle: str,
) -> CallableAlignment:
    reference = read_fasta(reference_path)
    intervals = read_bed(bed_path)
    sample_names = read_sample_names(sample_table, organelle)
    template, coordinate_to_index, sites = build_track_template(reference, intervals)

    raw_samples, raw_variants = read_vcf_variants(raw_vcf_path)
    filtered_samples, filtered_variants = read_vcf_variants(filtered_vcf_path)
    if raw_samples != sample_names:
        raise CallableConsensusError("Raw VCF sample order does not match included samples")
    if filtered_samples != sample_names:
        raise CallableConsensusError(
            "Filtered VCF sample order does not match included samples"
        )

    raw_site_keys = {
        (variant.chrom, variant.position)
        for variant in raw_variants
        if (variant.chrom, variant.position) in coordinate_to_index
    }
    filtered_site_keys = {
        (variant.chrom, variant.position)
        for variant in filtered_variants
        if (variant.chrom, variant.position) in coordinate_to_index
    }
    failed_site_indexes = {
        coordinate_to_index[site_key]
        for site_key in raw_site_keys - filtered_site_keys
    }

    sequence_parts = {sample: template.copy() for sample in sample_names}
    filtered_variant_sites = 0
    for variant in filtered_variants:
        index = coordinate_to_index.get((variant.chrom, variant.position))
        if index is None:
            continue
        if len(variant.ref) != 1 or len(variant.alt) != 1 or "," in variant.alt:
            continue
        filtered_variant_sites += 1
        for sample in sample_names:
            sequence_parts[sample][index] = genotype_to_base(
                variant.genotypes[sample],
                variant.ref,
                variant.alt,
            )

    for sample in sample_names:
        depth_path = depth_path_for_sample(depth_dir, sample)
        if not depth_path.exists():
            raise CallableConsensusError(f"Missing depth file for {sample}: {depth_path}")
        covered = read_covered_indexes(depth_path, coordinate_to_index, min_depth)
        sequence = sequence_parts[sample]
        for index in failed_site_indexes:
            sequence[index] = "N"
        for index in range(len(sequence)):
            if index not in covered:
                sequence[index] = "N"

    return CallableAlignment(
        sample_names=sample_names,
        sequences={sample: "".join(sequence_parts[sample]) for sample in sample_names},
        sites=sites,
        filtered_variant_sites=filtered_variant_sites,
        masked_failed_variant_sites=len(failed_site_indexes),
    )


def read_consensus_inputs(
    track_table: Path,
    variant_calling_dir: Path,
    variant_filtering_dir: Path,
    run_label: str,
) -> list[ConsensusInput]:
    tracks = {
        row["track_id"]: Path(row["bed_path"])
        for row in read_tsv(track_table)
        if row["track_id"] in {"cpdna_population_sites", "mtdna_high_confidence_unique"}
    }
    raw_by_organelle = {
        row["organelle"]: row
        for row in read_tsv(
            variant_calling_dir / labeled_output_name("variant_calling_summary.tsv", run_label)
        )
    }
    filtered_by_organelle = {
        row["organelle"]: row
        for row in read_tsv(
            variant_filtering_dir
            / labeled_output_name("variant_filtering_summary.tsv", run_label)
        )
    }
    inputs: list[ConsensusInput] = []
    for organelle, track_id in (
        ("cpDNA", "cpdna_population_sites"),
        ("mtDNA", "mtdna_high_confidence_unique"),
    ):
        if track_id not in tracks:
            raise CallableConsensusError(f"Missing track BED for {track_id}")
        if organelle not in raw_by_organelle:
            raise CallableConsensusError(f"Missing raw variant summary row for {organelle}")
        if organelle not in filtered_by_organelle:
            raise CallableConsensusError(
                f"Missing filtered variant summary row for {organelle}"
            )
        raw_row = raw_by_organelle[organelle]
        filtered_row = filtered_by_organelle[organelle]
        raw_vcf_path = Path(raw_row["raw_vcf_path"])
        filtered_vcf_path = Path(filtered_row["filtered_vcf_path"])
        if not raw_vcf_path.exists():
            raise CallableConsensusError(f"Missing raw VCF: {raw_vcf_path}")
        if not filtered_vcf_path.exists():
            raise CallableConsensusError(f"Missing filtered VCF: {filtered_vcf_path}")
        inputs.append(
            ConsensusInput(
                organelle=organelle,
                track_id=track_id,
                bed_path=tracks[track_id],
                raw_records=int(raw_row["variant_records"]),
                filtered_records=int(filtered_row["filtered_records"]),
                raw_vcf_path=raw_vcf_path,
                filtered_vcf_path=filtered_vcf_path,
            )
        )
    return inputs


def write_fasta(
    path: Path,
    sample_names: list[str],
    sequences: dict[str, str],
    line_width: int = 80,
) -> None:
    lines: list[str] = []
    for sample in sample_names:
        sequence = sequences[sample]
        lines.append(f">{sample}")
        for start in range(0, len(sequence), line_width):
            lines.append(sequence[start : start + line_width])
    path.write_text("\n".join(lines) + "\n")


def write_site_table(path: Path, sites: list[dict[str, str]]) -> None:
    write_tsv(path, sites, ["site_index", "chrom", "position", "reference_base"])


def consensus_output_paths(
    consensus_input: ConsensusInput,
    output_dir: Path,
    run_label: str,
) -> tuple[Path, Path]:
    label = f".{run_label}" if run_label else ""
    fasta = output_dir / f"{consensus_input.organelle}{label}.callable_consensus.fa"
    sites = output_dir / f"{consensus_input.organelle}{label}.callable_sites.tsv"
    return fasta, sites


def build_one_consensus(
    consensus_input: ConsensusInput,
    reference_path: Path,
    sample_table: Path,
    depth_dir: Path,
    output_dir: Path,
    run_label: str,
    min_depth: int,
) -> ConsensusResult:
    alignment = build_callable_consensus(
        reference_path=reference_path,
        bed_path=consensus_input.bed_path,
        sample_table=sample_table,
        depth_dir=depth_dir,
        raw_vcf_path=consensus_input.raw_vcf_path,
        filtered_vcf_path=consensus_input.filtered_vcf_path,
        min_depth=min_depth,
        organelle=consensus_input.organelle,
    )
    fasta_path, site_table_path = consensus_output_paths(
        consensus_input,
        output_dir,
        run_label,
    )
    write_fasta(fasta_path, alignment.sample_names, alignment.sequences)
    write_site_table(site_table_path, alignment.sites)
    return consensus_input.to_result(
        sample_count=len(alignment.sample_names),
        consensus_length=alignment.consensus_length,
        filtered_variant_sites=alignment.filtered_variant_sites,
        masked_failed_variant_sites=alignment.masked_failed_variant_sites,
        missing_bases=alignment.missing_bases,
        fasta_path=fasta_path,
        site_table_path=site_table_path,
    )


def run_callable_consensus(
    reference_path: Path = DEFAULT_REFERENCE,
    sample_table: Path = DEFAULT_SAMPLE_TABLE,
    track_table: Path = DEFAULT_TRACK_TABLE,
    depth_dir: Path = DEFAULT_DEPTH_DIR,
    variant_calling_dir: Path = DEFAULT_VARIANT_CALLING_DIR,
    variant_filtering_dir: Path = DEFAULT_VARIANT_FILTERING_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
    min_depth: int = DEFAULT_MIN_DEPTH,
) -> list[ConsensusResult]:
    output_dir.mkdir(parents=True, exist_ok=True)
    inputs = read_consensus_inputs(
        track_table=track_table,
        variant_calling_dir=variant_calling_dir,
        variant_filtering_dir=variant_filtering_dir,
        run_label=run_label,
    )
    results = [
        build_one_consensus(
            consensus_input=consensus_input,
            reference_path=reference_path,
            sample_table=sample_table,
            depth_dir=depth_dir,
            output_dir=output_dir,
            run_label=run_label,
            min_depth=min_depth,
        )
        for consensus_input in inputs
    ]
    write_consensus_outputs(output_dir, results, run_label=run_label, min_depth=min_depth)
    return results


def write_consensus_outputs(
    output_dir: Path,
    results: list[ConsensusResult],
    run_label: str,
    min_depth: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / labeled_output_name("callable_consensus_summary.tsv", run_label),
        [
            {
                "organelle": result.organelle,
                "track_id": result.track_id,
                "sample_count": str(result.sample_count),
                "consensus_length": str(result.consensus_length),
                "raw_records": str(result.raw_records),
                "filtered_records": str(result.filtered_records),
                "filtered_variant_sites": str(result.filtered_variant_sites),
                "masked_failed_variant_sites": str(result.masked_failed_variant_sites),
                "missing_bases": str(result.missing_bases),
                "bed_path": result.bed_path.as_posix(),
                "raw_vcf_path": result.raw_vcf_path.as_posix(),
                "filtered_vcf_path": result.filtered_vcf_path.as_posix(),
                "alignment_fasta_path": result.fasta_path.as_posix(),
                "site_table_path": result.site_table_path.as_posix(),
            }
            for result in results
        ],
        [
            "organelle",
            "track_id",
            "sample_count",
            "consensus_length",
            "raw_records",
            "filtered_records",
            "filtered_variant_sites",
            "masked_failed_variant_sites",
            "missing_bases",
            "bed_path",
            "raw_vcf_path",
            "filtered_vcf_path",
            "alignment_fasta_path",
            "site_table_path",
        ],
    )
    write_report(
        output_dir / labeled_output_name("callable_consensus_report.md", run_label),
        results=results,
        run_label=run_label,
        min_depth=min_depth,
    )


def write_report(
    path: Path,
    results: list[ConsensusResult],
    run_label: str,
    min_depth: int,
) -> None:
    label = run_label or "full"
    lines = [
        "# Step 10 Callable-Site Consensus Alignment",
        "",
        "This step builds full callable-site FASTA alignments for cpDNA and mtDNA.",
        "Each alignment follows the Step 4 population-genetic BED track, starts",
        "from the annotated organelle reference, overlays Step 8 filtered haploid",
        "SNP genotypes, masks Step 7 raw variant sites that failed filtering, and",
        "uses Step 5 depth files to write `N` at bases below the minimum depth.",
        "",
        "## Run",
        "",
        f"- Run label: `{label}`",
        f"- Minimum depth for a non-missing consensus base: `{min_depth}`",
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
                f"- Consensus length: {result.consensus_length}",
                f"- Raw variant records considered: {result.raw_records}",
                f"- Filtered SNP records available: {result.filtered_records}",
                f"- Filtered SNP sites applied inside track: {result.filtered_variant_sites}",
                f"- Raw-only failed variant sites masked: {result.masked_failed_variant_sites}",
                f"- Missing consensus bases: {result.missing_bases}",
                f"- FASTA: `{result.fasta_path}`",
                f"- Site table: `{result.site_table_path}`",
                "",
            ]
        )
    path.write_text("\n".join(lines))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Step 10: build full callable-site cpDNA/mtDNA consensus alignments."
    )
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--sample-table", type=Path, default=DEFAULT_SAMPLE_TABLE)
    parser.add_argument("--track-table", type=Path, default=DEFAULT_TRACK_TABLE)
    parser.add_argument("--depth-dir", type=Path, default=DEFAULT_DEPTH_DIR)
    parser.add_argument("--variant-calling-dir", type=Path, default=DEFAULT_VARIANT_CALLING_DIR)
    parser.add_argument(
        "--variant-filtering-dir",
        type=Path,
        default=DEFAULT_VARIANT_FILTERING_DIR,
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    parser.add_argument("--min-depth", type=int, default=DEFAULT_MIN_DEPTH)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    results = run_callable_consensus(
        reference_path=args.reference,
        sample_table=args.sample_table,
        track_table=args.track_table,
        depth_dir=args.depth_dir,
        variant_calling_dir=args.variant_calling_dir,
        variant_filtering_dir=args.variant_filtering_dir,
        output_dir=args.output_dir,
        run_label=args.run_label,
        min_depth=args.min_depth,
    )
    for result in results:
        print(
            f"{result.organelle}: {result.consensus_length} callable consensus sites "
            f"across {result.sample_count} samples"
        )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
