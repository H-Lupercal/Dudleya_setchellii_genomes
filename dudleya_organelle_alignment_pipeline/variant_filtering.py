"""Filter raw haploid cpDNA/mtDNA variants for downstream analyses.

This is Step 8 of the pipeline. It consumes the Step 7 raw VCF summary and
writes filtered biallelic SNP VCFs for PCA, trees, Fst, and clustering inputs.
Consensus FASTA generation happens in a later step.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from dudleya_organelle_alignment_pipeline.pilot_alignment import shlex_join
from dudleya_organelle_alignment_pipeline.variant_calling import (
    count_vcf_records,
    labeled_output_name,
)


DEFAULT_VARIANT_CALLING_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/08_variant_calling"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/09_variant_filtering"
)
DEFAULT_RUN_LABEL = "primary"
DEFAULT_MAX_MISSING_FRACTION = 0.2
DEFAULT_MIN_MINOR_ALLELE_COUNT = 2
DEFAULT_THREADS = 4


class VariantFilteringError(RuntimeError):
    """Raised when Step 8 cannot safely filter variants."""


@dataclass(frozen=True)
class FilterInput:
    organelle: str
    track_id: str
    sample_count: int
    raw_records: int
    raw_vcf_path: Path
    raw_vcf_index_path: Path

    @property
    def filtered_vcf_path(self) -> Path:
        return Path(str(self.raw_vcf_path).replace(".raw.vcf.gz", ".filtered.vcf.gz"))

    @property
    def filtered_index_path(self) -> Path:
        return Path(f"{self.filtered_vcf_path}.tbi")

    @property
    def log_path(self) -> Path:
        return Path(str(self.raw_vcf_path).replace(".raw.vcf.gz", ".filtered.bcftools.log"))


@dataclass(frozen=True)
class FilterResult:
    organelle: str
    track_id: str
    sample_count: int
    raw_records: int
    filtered_records: int
    raw_vcf_path: Path
    filtered_vcf_path: Path
    filtered_index_path: Path
    log_path: Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def require_bcftools() -> None:
    if shutil.which("bcftools") is None:
        raise VariantFilteringError(
            "Missing required tool: bcftools. Activate the pipeline environment first."
        )


def read_filter_inputs(summary_path: Path) -> list[FilterInput]:
    inputs: list[FilterInput] = []
    for row in read_tsv(summary_path):
        raw_vcf_path = Path(row["raw_vcf_path"])
        raw_vcf_index_path = Path(row["raw_vcf_index_path"])
        if not raw_vcf_path.exists():
            raise VariantFilteringError(f"Missing raw VCF: {raw_vcf_path}")
        if not raw_vcf_index_path.exists():
            raise VariantFilteringError(f"Missing raw VCF index: {raw_vcf_index_path}")
        inputs.append(
            FilterInput(
                organelle=row["organelle"],
                track_id=row["track_id"],
                sample_count=int(row["sample_count"]),
                raw_records=int(row["variant_records"]),
                raw_vcf_path=raw_vcf_path,
                raw_vcf_index_path=raw_vcf_index_path,
            )
        )
    if not inputs:
        raise VariantFilteringError(f"No variant-call rows found in {summary_path}")
    return inputs


def output_filtered_vcf_path(filter_input: FilterInput, output_dir: Path) -> Path:
    return output_dir / filter_input.filtered_vcf_path.name


def output_log_path(filter_input: FilterInput, output_dir: Path) -> Path:
    return output_dir / filter_input.log_path.name


def build_filter_command(
    raw_vcf: Path,
    filtered_vcf: Path,
    max_missing_fraction: float,
    min_minor_allele_count: int,
    threads: int,
) -> list[str]:
    return [
        "bcftools",
        "view",
        "--threads",
        str(threads),
        "-m2",
        "-M2",
        "-v",
        "snps",
        "--min-ac",
        f"{min_minor_allele_count}:minor",
        "-i",
        f"F_MISSING<={max_missing_fraction:g}",
        "-Oz",
        "-o",
        filtered_vcf.as_posix(),
        raw_vcf.as_posix(),
    ]


def build_index_command(filtered_vcf: Path) -> list[str]:
    return ["bcftools", "index", "-t", filtered_vcf.as_posix()]


def outputs_are_ready(filtered_vcf: Path, filtered_index: Path) -> bool:
    return (
        filtered_vcf.exists()
        and filtered_vcf.stat().st_size > 0
        and filtered_index.exists()
        and filtered_index.stat().st_size > 0
    )


def filter_one_input(
    filter_input: FilterInput,
    output_dir: Path,
    max_missing_fraction: float,
    min_minor_allele_count: int,
    threads: int,
    force: bool = False,
) -> tuple[FilterResult, list[dict[str, str]]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    filtered_vcf = output_filtered_vcf_path(filter_input, output_dir)
    filtered_index = Path(f"{filtered_vcf}.tbi")
    log_path = output_log_path(filter_input, output_dir)
    filter_command = build_filter_command(
        raw_vcf=filter_input.raw_vcf_path,
        filtered_vcf=filtered_vcf,
        max_missing_fraction=max_missing_fraction,
        min_minor_allele_count=min_minor_allele_count,
        threads=threads,
    )
    index_command = build_index_command(filtered_vcf)
    command_rows = [
        {
            "organelle": filter_input.organelle,
            "track_id": filter_input.track_id,
            "step": "filter",
            "command": shlex_join(filter_command),
        },
        {
            "organelle": filter_input.organelle,
            "track_id": filter_input.track_id,
            "step": "index",
            "command": shlex_join(index_command),
        },
    ]
    if outputs_are_ready(filtered_vcf, filtered_index) and not force:
        command_rows.append(
            {
                "organelle": filter_input.organelle,
                "track_id": filter_input.track_id,
                "step": "reuse_existing_outputs",
                "command": "outputs already present; pass --force to regenerate",
            }
        )
        return (
            FilterResult(
                organelle=filter_input.organelle,
                track_id=filter_input.track_id,
                sample_count=filter_input.sample_count,
                raw_records=filter_input.raw_records,
                filtered_records=count_vcf_records(filtered_vcf),
                raw_vcf_path=filter_input.raw_vcf_path,
                filtered_vcf_path=filtered_vcf,
                filtered_index_path=filtered_index,
                log_path=log_path,
            ),
            command_rows,
        )

    with log_path.open("w") as log_handle:
        for command, step_name in (
            (filter_command, "filter"),
            (index_command, "index"),
        ):
            completed = subprocess.run(
                command,
                stdout=subprocess.DEVNULL,
                stderr=log_handle,
                check=False,
            )
            if completed.returncode:
                raise VariantFilteringError(
                    f"bcftools {step_name} failed for {filter_input.organelle}; "
                    f"see {log_path}"
                )

    return (
        FilterResult(
            organelle=filter_input.organelle,
            track_id=filter_input.track_id,
            sample_count=filter_input.sample_count,
            raw_records=filter_input.raw_records,
            filtered_records=count_vcf_records(filtered_vcf),
            raw_vcf_path=filter_input.raw_vcf_path,
            filtered_vcf_path=filtered_vcf,
            filtered_index_path=filtered_index,
            log_path=log_path,
        ),
        command_rows,
    )


def run_variant_filtering(
    variant_calling_dir: Path = DEFAULT_VARIANT_CALLING_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
    max_missing_fraction: float = DEFAULT_MAX_MISSING_FRACTION,
    min_minor_allele_count: int = DEFAULT_MIN_MINOR_ALLELE_COUNT,
    threads: int = DEFAULT_THREADS,
    force: bool = False,
) -> list[FilterResult]:
    require_bcftools()
    summary_path = variant_calling_dir / labeled_output_name(
        "variant_calling_summary.tsv", run_label
    )
    filter_inputs = read_filter_inputs(summary_path)
    results: list[FilterResult] = []
    command_rows: list[dict[str, str]] = []
    for filter_input in filter_inputs:
        result, rows = filter_one_input(
            filter_input=filter_input,
            output_dir=output_dir,
            max_missing_fraction=max_missing_fraction,
            min_minor_allele_count=min_minor_allele_count,
            threads=threads,
            force=force,
        )
        results.append(result)
        command_rows.extend(rows)
    write_filtering_outputs(
        output_dir=output_dir,
        results=results,
        command_rows=command_rows,
        run_label=run_label,
        max_missing_fraction=max_missing_fraction,
        min_minor_allele_count=min_minor_allele_count,
    )
    return results


def write_filtering_outputs(
    output_dir: Path,
    results: list[FilterResult],
    command_rows: list[dict[str, str]],
    run_label: str,
    max_missing_fraction: float,
    min_minor_allele_count: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / labeled_output_name("filtering.commands.tsv", run_label),
        command_rows,
        ["organelle", "track_id", "step", "command"],
    )
    write_tsv(
        output_dir / labeled_output_name("variant_filtering_summary.tsv", run_label),
        [
            {
                "organelle": result.organelle,
                "track_id": result.track_id,
                "sample_count": str(result.sample_count),
                "raw_records": str(result.raw_records),
                "filtered_records": str(result.filtered_records),
                "raw_vcf_path": result.raw_vcf_path.as_posix(),
                "filtered_vcf_path": result.filtered_vcf_path.as_posix(),
                "filtered_vcf_index_path": result.filtered_index_path.as_posix(),
                "log_path": result.log_path.as_posix(),
            }
            for result in results
        ],
        [
            "organelle",
            "track_id",
            "sample_count",
            "raw_records",
            "filtered_records",
            "raw_vcf_path",
            "filtered_vcf_path",
            "filtered_vcf_index_path",
            "log_path",
        ],
    )
    write_report(
        output_dir / labeled_output_name("variant_filtering_report.md", run_label),
        results=results,
        run_label=run_label,
        max_missing_fraction=max_missing_fraction,
        min_minor_allele_count=min_minor_allele_count,
    )


def write_report(
    path: Path,
    results: list[FilterResult],
    run_label: str,
    max_missing_fraction: float,
    min_minor_allele_count: int,
) -> None:
    label = run_label or "full"
    lines = [
        "# Step 8 Variant Filtering",
        "",
        "This step filters the raw haploid cpDNA and mtDNA variant calls from Step 7.",
        "Consensus FASTA generation, alignments, PCA, and trees happen in later steps.",
        "",
        "## Run",
        "",
        f"- Run label: `{label}`",
        "- Variant type retained: biallelic SNPs",
        f"- Maximum missing genotype fraction: {max_missing_fraction:g}",
        f"- Minimum minor allele count: {min_minor_allele_count}",
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
                f"- Raw records: {result.raw_records}",
                f"- Filtered records: {result.filtered_records}",
                f"- Filtered VCF: `{result.filtered_vcf_path}`",
                f"- Index: `{result.filtered_index_path}`",
                f"- Log: `{result.log_path}`",
                "",
            ]
        )
    path.write_text("\n".join(lines))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Step 8: filter raw haploid cpDNA/mtDNA variants."
    )
    parser.add_argument("--variant-calling-dir", type=Path, default=DEFAULT_VARIANT_CALLING_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    parser.add_argument(
        "--max-missing-fraction",
        type=float,
        default=DEFAULT_MAX_MISSING_FRACTION,
    )
    parser.add_argument(
        "--min-minor-allele-count",
        type=int,
        default=DEFAULT_MIN_MINOR_ALLELE_COUNT,
    )
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    results = run_variant_filtering(
        variant_calling_dir=args.variant_calling_dir,
        output_dir=args.output_dir,
        run_label=args.run_label,
        max_missing_fraction=args.max_missing_fraction,
        min_minor_allele_count=args.min_minor_allele_count,
        threads=args.threads,
        force=args.force,
    )
    for result in results:
        print(
            f"{result.organelle}: {result.filtered_records} filtered records "
            f"from {result.raw_records} raw records"
        )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
