"""Build SNP-only cpDNA/mtDNA FASTA alignments from filtered haploid VCFs.

This stage converts the filtered biallelic SNP VCFs
into all-sample FASTA alignments suitable for quick phylogenetic trees and PCA
matrix generation. Reference-length consensus FASTAs can be generated later if
needed for full callable-site analyses.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from dudleya_organelle_alignment_pipeline.variant_calling import labeled_output_name


DEFAULT_VARIANT_FILTERING_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/09_variant_filtering"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/10_snp_alignment"
)
DEFAULT_RUN_LABEL = "primary"


class SnpAlignmentError(RuntimeError):
    """Raised when this stage cannot safely build SNP alignments."""


@dataclass(frozen=True)
class SnpAlignment:
    sample_names: list[str]
    sequences: dict[str, str]
    sites: list[dict[str, str]]

    @property
    def site_count(self) -> int:
        return len(self.sites)

    @property
    def missing_bases(self) -> int:
        return sum(sequence.count("N") for sequence in self.sequences.values())


@dataclass(frozen=True)
class SnpAlignmentInput:
    organelle: str
    track_id: str
    sample_count: int
    filtered_records: int
    filtered_vcf_path: Path
    filtered_vcf_index_path: Path

    def to_result(
        self,
        alignment: SnpAlignment,
        alignment_fasta_path: Path,
        site_table_path: Path,
    ) -> "SnpAlignmentResult":
        return SnpAlignmentResult(
            organelle=self.organelle,
            track_id=self.track_id,
            sample_count=len(alignment.sample_names),
            filtered_records=self.filtered_records,
            alignment_sites=alignment.site_count,
            missing_bases=alignment.missing_bases,
            filtered_vcf_path=self.filtered_vcf_path,
            alignment_fasta_path=alignment_fasta_path,
            site_table_path=site_table_path,
        )


@dataclass(frozen=True)
class SnpAlignmentResult:
    organelle: str
    track_id: str
    sample_count: int
    filtered_records: int
    alignment_sites: int
    missing_bases: int
    filtered_vcf_path: Path
    alignment_fasta_path: Path
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


def read_alignment_inputs(summary_path: Path) -> list[SnpAlignmentInput]:
    inputs: list[SnpAlignmentInput] = []
    for row in read_tsv(summary_path):
        filtered_vcf_path = Path(row["filtered_vcf_path"])
        filtered_vcf_index_path = Path(row["filtered_vcf_index_path"])
        if not filtered_vcf_path.exists():
            raise SnpAlignmentError(f"Missing filtered VCF: {filtered_vcf_path}")
        if not filtered_vcf_index_path.exists():
            raise SnpAlignmentError(f"Missing filtered VCF index: {filtered_vcf_index_path}")
        inputs.append(
            SnpAlignmentInput(
                organelle=row["organelle"],
                track_id=row["track_id"],
                sample_count=int(row["sample_count"]),
                filtered_records=int(row["filtered_records"]),
                filtered_vcf_path=filtered_vcf_path,
                filtered_vcf_index_path=filtered_vcf_index_path,
            )
        )
    if not inputs:
        raise SnpAlignmentError(f"No filtered variant rows found in {summary_path}")
    return inputs


def genotype_to_base(genotype_field: str, ref: str, alt: str) -> str:
    genotype = genotype_field.split(":", 1)[0]
    allele = genotype.replace("|", "/").split("/", 1)[0]
    if allele == "0":
        return ref.upper()
    if allele == "1":
        return alt.upper()
    return "N"


def build_snp_alignment(vcf_path: Path) -> SnpAlignment:
    sample_names: list[str] | None = None
    sequence_parts: dict[str, list[str]] = {}
    sites: list[dict[str, str]] = []
    with open_text(vcf_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                columns = line.split("\t")
                sample_names = columns[9:]
                sequence_parts = {sample: [] for sample in sample_names}
                continue
            if sample_names is None:
                raise SnpAlignmentError(f"VCF header missing #CHROM line: {vcf_path}")
            columns = line.split("\t")
            chrom, position, ref, alt = columns[0], columns[1], columns[3], columns[4]
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                raise SnpAlignmentError(
                    f"Expected biallelic SNP in filtered VCF, found {chrom}:{position}"
                )
            sites.append(
                {
                    "site_index": str(len(sites) + 1),
                    "chrom": chrom,
                    "position": position,
                    "ref": ref,
                    "alt": alt,
                }
            )
            for sample, genotype_field in zip(sample_names, columns[9:], strict=True):
                sequence_parts[sample].append(genotype_to_base(genotype_field, ref, alt))

    if sample_names is None:
        raise SnpAlignmentError(f"VCF header missing #CHROM line: {vcf_path}")
    return SnpAlignment(
        sample_names=sample_names,
        sequences={sample: "".join(sequence_parts[sample]) for sample in sample_names},
        sites=sites,
    )


def write_fasta(path: Path, alignment: SnpAlignment, line_width: int = 80) -> None:
    lines: list[str] = []
    for sample in alignment.sample_names:
        sequence = alignment.sequences[sample]
        lines.append(f">{sample}")
        for start in range(0, len(sequence), line_width):
            lines.append(sequence[start : start + line_width])
    path.write_text("\n".join(lines) + "\n")


def write_site_table(path: Path, alignment: SnpAlignment) -> None:
    write_tsv(path, alignment.sites, ["site_index", "chrom", "position", "ref", "alt"])


def alignment_output_paths(
    alignment_input: SnpAlignmentInput,
    output_dir: Path,
    run_label: str,
) -> tuple[Path, Path]:
    label = f".{run_label}" if run_label else ""
    fasta = output_dir / f"{alignment_input.organelle}{label}.snp_alignment.fa"
    sites = output_dir / f"{alignment_input.organelle}{label}.snp_sites.tsv"
    return fasta, sites


def build_one_alignment(
    alignment_input: SnpAlignmentInput,
    output_dir: Path,
    run_label: str,
) -> SnpAlignmentResult:
    output_dir.mkdir(parents=True, exist_ok=True)
    alignment = build_snp_alignment(alignment_input.filtered_vcf_path)
    fasta_path, site_table_path = alignment_output_paths(alignment_input, output_dir, run_label)
    write_fasta(fasta_path, alignment)
    write_site_table(site_table_path, alignment)
    return alignment_input.to_result(
        alignment=alignment,
        alignment_fasta_path=fasta_path,
        site_table_path=site_table_path,
    )


def run_snp_alignment(
    variant_filtering_dir: Path = DEFAULT_VARIANT_FILTERING_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
) -> list[SnpAlignmentResult]:
    summary_path = variant_filtering_dir / labeled_output_name(
        "variant_filtering_summary.tsv", run_label
    )
    inputs = read_alignment_inputs(summary_path)
    results = [
        build_one_alignment(
            alignment_input=alignment_input,
            output_dir=output_dir,
            run_label=run_label,
        )
        for alignment_input in inputs
    ]
    write_alignment_outputs(output_dir=output_dir, results=results, run_label=run_label)
    return results


def write_alignment_outputs(
    output_dir: Path,
    results: list[SnpAlignmentResult],
    run_label: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / labeled_output_name("snp_alignment_summary.tsv", run_label),
        [
            {
                "organelle": result.organelle,
                "track_id": result.track_id,
                "sample_count": str(result.sample_count),
                "filtered_records": str(result.filtered_records),
                "alignment_sites": str(result.alignment_sites),
                "missing_bases": str(result.missing_bases),
                "filtered_vcf_path": result.filtered_vcf_path.as_posix(),
                "alignment_fasta_path": result.alignment_fasta_path.as_posix(),
                "site_table_path": result.site_table_path.as_posix(),
            }
            for result in results
        ],
        [
            "organelle",
            "track_id",
            "sample_count",
            "filtered_records",
            "alignment_sites",
            "missing_bases",
            "filtered_vcf_path",
            "alignment_fasta_path",
            "site_table_path",
        ],
    )
    write_report(
        output_dir / labeled_output_name("snp_alignment_report.md", run_label),
        results=results,
        run_label=run_label,
    )


def write_report(path: Path, results: list[SnpAlignmentResult], run_label: str) -> None:
    label = run_label or "full"
    lines = [
        "# SNP Alignment",
        "",
        "This step converts filtered haploid cpDNA and mtDNA SNP VCFs into",
        "SNP-only FASTA alignments. These alignments are intended for quick",
        "tree-building and matrix-based analyses; full reference-length",
        "consensus FASTAs can be generated in a later step if needed.",
        "",
        "## Run",
        "",
        f"- Run label: `{label}`",
        "- Haploid genotype encoding: `0` uses REF, `1` uses ALT, missing uses `N`.",
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
                f"- Filtered records: {result.filtered_records}",
                f"- Alignment sites: {result.alignment_sites}",
                f"- Missing alignment bases: {result.missing_bases}",
                f"- FASTA: `{result.alignment_fasta_path}`",
                f"- Site table: `{result.site_table_path}`",
                "",
            ]
        )
    path.write_text("\n".join(lines))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build SNP-only cpDNA/mtDNA FASTA alignments."
    )
    parser.add_argument("--variant-filtering-dir", type=Path, default=DEFAULT_VARIANT_FILTERING_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    results = run_snp_alignment(
        variant_filtering_dir=args.variant_filtering_dir,
        output_dir=args.output_dir,
        run_label=args.run_label,
    )
    for result in results:
        print(
            f"{result.organelle}: {result.alignment_sites} SNP-alignment sites "
            f"across {result.sample_count} samples"
        )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
