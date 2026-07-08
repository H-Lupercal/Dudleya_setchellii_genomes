"""Run haploid cpDNA/mtDNA variant calling for downstream analyses.

This stage uses the downstream sample set and
the population-genetic tracks to call raw haploid variants separately
for cpDNA and mtDNA. Filtering and consensus generation happen in later steps.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from dudleya_organelle_alignment_pipeline.pilot_alignment import safe_sample_name, shlex_join


DEFAULT_SAMPLE_TABLE = Path(
    "dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/"
    "included_samples.tsv"
)
DEFAULT_TRACK_TABLE = Path(
    "dudleya_organelle_alignment_pipeline/results/05_analysis_masks/analysis_tracks.tsv"
)
DEFAULT_REFERENCE = Path(
    "dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa"
)
DEFAULT_BAM_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/bam"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/08_variant_calling"
)
DEFAULT_MIN_MAPQ = 20
DEFAULT_MIN_BASEQ = 20
DEFAULT_MAX_DEPTH = 10000
DEFAULT_THREADS = 4


class VariantCallingError(RuntimeError):
    """Raised when this stage cannot safely call variants."""


@dataclass(frozen=True)
class VariantSample:
    sample_id: str
    safe_sample_id: str
    row: dict[str, str]
    bam_path: Path
    bam_index_path: Path


@dataclass(frozen=True)
class VariantTrack:
    organelle: str
    track_id: str
    bed_path: Path
    output_prefix: str


@dataclass(frozen=True)
class VariantCallInputs:
    bam_list_path: Path
    sample_names_path: Path
    final_vcf_path: Path
    final_index_path: Path
    pre_reheader_vcf_path: Path
    log_path: Path


@dataclass(frozen=True)
class VariantCallResult:
    organelle: str
    track_id: str
    sample_count: int
    variant_records: int
    final_vcf_path: Path
    final_index_path: Path
    log_path: Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_variant_samples(
    sample_table: Path,
    bam_dir: Path,
    sample_limit: int | None = None,
    sample_ids: set[str] | None = None,
) -> list[VariantSample]:
    samples: list[VariantSample] = []
    for row in read_tsv(sample_table):
        sample_id = row["sample_id"]
        if sample_ids is not None and sample_id not in sample_ids:
            continue
        if row.get("downstream_cpDNA_use") != "include":
            continue
        if row.get("downstream_mtDNA_use") != "include":
            continue
        safe_id = safe_sample_name(sample_id)
        bam_path = bam_dir / f"{safe_id}.organelle.sorted.bam"
        bam_index_path = bam_dir / f"{safe_id}.organelle.sorted.bam.bai"
        if not bam_path.exists():
            raise VariantCallingError(f"Missing BAM for {sample_id}: {bam_path}")
        if not bam_index_path.exists():
            raise VariantCallingError(f"Missing BAM index for {sample_id}: {bam_index_path}")
        samples.append(
            VariantSample(
                sample_id=sample_id,
                safe_sample_id=safe_id,
                row=row,
                bam_path=bam_path,
                bam_index_path=bam_index_path,
            )
        )
        if sample_limit is not None and len(samples) >= sample_limit:
            break

    if sample_ids is not None:
        found = {sample.sample_id for sample in samples}
        missing = sorted(sample_ids - found)
        if missing:
            raise VariantCallingError(
                "Requested samples are absent from the downstream included set: "
                + ", ".join(missing)
            )
    if not samples:
        raise VariantCallingError(f"No downstream variant samples found in {sample_table}")
    return samples


def read_variant_tracks(track_table: Path, run_label: str = "") -> dict[str, VariantTrack]:
    tracks: dict[str, VariantTrack] = {}
    output_prefixes = {
        "cpDNA": label_output_prefix("cpDNA.raw", run_label),
        "mtDNA": label_output_prefix("mtDNA.raw", run_label),
    }
    wanted = {
        "cpdna_population_sites": "cpDNA",
        "mtdna_high_confidence_unique": "mtDNA",
    }
    for row in read_tsv(track_table):
        track_id = row["track_id"]
        if track_id not in wanted:
            continue
        organelle = wanted[track_id]
        bed_path = Path(row["bed_path"])
        if not bed_path.exists():
            raise VariantCallingError(f"Missing BED for {track_id}: {bed_path}")
        tracks[organelle] = VariantTrack(
            organelle=organelle,
            track_id=track_id,
            bed_path=bed_path,
            output_prefix=output_prefixes[organelle],
        )

    missing = sorted(set(output_prefixes) - set(tracks))
    if missing:
        raise VariantCallingError(
            "Missing population-genetic variant tracks for: " + ", ".join(missing)
        )
    return tracks


def label_output_prefix(base_prefix: str, run_label: str) -> str:
    if not run_label:
        return base_prefix
    return f"{base_prefix.rsplit('.', 1)[0]}.{run_label}.{base_prefix.rsplit('.', 1)[1]}"


def write_variant_call_inputs(
    track: VariantTrack,
    samples: list[VariantSample],
    output_dir: Path,
) -> VariantCallInputs:
    output_dir.mkdir(parents=True, exist_ok=True)
    bam_list_path = output_dir / f"{track.output_prefix}.bam_list.txt"
    sample_names_path = output_dir / f"{track.output_prefix}.sample_names.txt"
    final_vcf_path = output_dir / f"{track.output_prefix}.vcf.gz"
    pre_reheader_vcf_path = output_dir / f"{track.output_prefix}.pre_reheader.vcf.gz"
    log_path = output_dir / f"{track.output_prefix}.bcftools.log"
    bam_list_path.write_text(
        "".join(f"{sample.bam_path.as_posix()}\n" for sample in samples)
    )
    sample_names_path.write_text("".join(f"{sample.sample_id}\n" for sample in samples))
    return VariantCallInputs(
        bam_list_path=bam_list_path,
        sample_names_path=sample_names_path,
        final_vcf_path=final_vcf_path,
        final_index_path=Path(f"{final_vcf_path}.tbi"),
        pre_reheader_vcf_path=pre_reheader_vcf_path,
        log_path=log_path,
    )


def build_bcftools_commands(
    track: VariantTrack,
    inputs: VariantCallInputs,
    reference: Path,
    min_mapq: int,
    min_baseq: int,
    max_depth: int,
    threads: int,
) -> dict[str, list[str]]:
    return {
        "mpileup": [
            "bcftools",
            "mpileup",
            "-Ou",
            "--threads",
            str(threads),
            "--ignore-RG",
            "--max-depth",
            str(max_depth),
            "-q",
            str(min_mapq),
            "-Q",
            str(min_baseq),
            "-a",
            "FORMAT/DP,FORMAT/AD",
            "-f",
            reference.as_posix(),
            "-R",
            track.bed_path.as_posix(),
            "-b",
            inputs.bam_list_path.as_posix(),
        ],
        "call": [
            "bcftools",
            "call",
            "--threads",
            str(threads),
            "--ploidy",
            "1",
            "-m",
            "-v",
            "-Oz",
            "-o",
            inputs.pre_reheader_vcf_path.as_posix(),
        ],
        "reheader": [
            "bcftools",
            "reheader",
            "-N",
            inputs.sample_names_path.as_posix(),
            "-o",
            inputs.final_vcf_path.as_posix(),
            inputs.pre_reheader_vcf_path.as_posix(),
        ],
        "index": ["bcftools", "index", "-t", inputs.final_vcf_path.as_posix()],
    }


def require_bcftools() -> None:
    if shutil.which("bcftools") is None:
        raise VariantCallingError(
            "Missing required tool: bcftools. Activate the pipeline environment first."
        )


def run_variant_call_for_track(
    track: VariantTrack,
    samples: list[VariantSample],
    output_dir: Path,
    reference: Path,
    min_mapq: int,
    min_baseq: int,
    max_depth: int,
    threads: int,
    force: bool = False,
) -> tuple[VariantCallResult, list[dict[str, str]]]:
    inputs = write_variant_call_inputs(track, samples, output_dir)
    commands = build_bcftools_commands(
        track=track,
        inputs=inputs,
        reference=reference,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        max_depth=max_depth,
        threads=threads,
    )
    command_rows = [
        {
            "organelle": track.organelle,
            "track_id": track.track_id,
            "step": "mpileup_call",
            "command": f"{shlex_join(commands['mpileup'])} | {shlex_join(commands['call'])}",
        },
        {
            "organelle": track.organelle,
            "track_id": track.track_id,
            "step": "reheader",
            "command": shlex_join(commands["reheader"]),
        },
        {
            "organelle": track.organelle,
            "track_id": track.track_id,
            "step": "index",
            "command": shlex_join(commands["index"]),
        },
    ]
    if outputs_are_ready(inputs) and not force:
        result = VariantCallResult(
            organelle=track.organelle,
            track_id=track.track_id,
            sample_count=len(samples),
            variant_records=count_vcf_records(inputs.final_vcf_path),
            final_vcf_path=inputs.final_vcf_path,
            final_index_path=inputs.final_index_path,
            log_path=inputs.log_path,
        )
        command_rows.append(
            {
                "organelle": track.organelle,
                "track_id": track.track_id,
                "step": "reuse_existing_outputs",
                "command": "outputs already present; pass --force to regenerate",
            }
        )
        return result, command_rows

    with inputs.log_path.open("w") as log_handle:
        mpileup_proc = subprocess.Popen(
            commands["mpileup"],
            stdout=subprocess.PIPE,
            stderr=log_handle,
        )
        assert mpileup_proc.stdout is not None
        call_proc = subprocess.Popen(
            commands["call"],
            stdin=mpileup_proc.stdout,
            stdout=subprocess.DEVNULL,
            stderr=log_handle,
        )
        mpileup_proc.stdout.close()
        call_return = call_proc.wait()
        mpileup_return = mpileup_proc.wait()
        if mpileup_return or call_return:
            raise VariantCallingError(
                f"bcftools mpileup/call failed for {track.organelle}; see {inputs.log_path}"
            )
        for step in ("reheader", "index"):
            completed = subprocess.run(
                commands[step],
                stdout=subprocess.DEVNULL,
                stderr=log_handle,
                check=False,
            )
            if completed.returncode:
                raise VariantCallingError(
                    f"bcftools {step} failed for {track.organelle}; see {inputs.log_path}"
                )

    inputs.pre_reheader_vcf_path.unlink(missing_ok=True)
    result = VariantCallResult(
        organelle=track.organelle,
        track_id=track.track_id,
        sample_count=len(samples),
        variant_records=count_vcf_records(inputs.final_vcf_path),
        final_vcf_path=inputs.final_vcf_path,
        final_index_path=inputs.final_index_path,
        log_path=inputs.log_path,
    )
    return result, command_rows


def outputs_are_ready(inputs: VariantCallInputs) -> bool:
    return (
        inputs.final_vcf_path.exists()
        and inputs.final_vcf_path.stat().st_size > 0
        and inputs.final_index_path.exists()
        and inputs.final_index_path.stat().st_size > 0
    )


def count_vcf_records(vcf_path: Path) -> int:
    count = 0
    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                count += 1
    return count


def run_variant_calling(
    sample_table: Path = DEFAULT_SAMPLE_TABLE,
    track_table: Path = DEFAULT_TRACK_TABLE,
    reference: Path = DEFAULT_REFERENCE,
    bam_dir: Path = DEFAULT_BAM_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    sample_limit: int | None = None,
    sample_ids: set[str] | None = None,
    run_label: str = "",
    min_mapq: int = DEFAULT_MIN_MAPQ,
    min_baseq: int = DEFAULT_MIN_BASEQ,
    max_depth: int = DEFAULT_MAX_DEPTH,
    threads: int = DEFAULT_THREADS,
    force: bool = False,
) -> list[VariantCallResult]:
    require_bcftools()
    samples = read_variant_samples(
        sample_table=sample_table,
        bam_dir=bam_dir,
        sample_limit=sample_limit,
        sample_ids=sample_ids,
    )
    tracks = read_variant_tracks(track_table, run_label=run_label)
    all_command_rows: list[dict[str, str]] = []
    results: list[VariantCallResult] = []
    for organelle in ("cpDNA", "mtDNA"):
        result, command_rows = run_variant_call_for_track(
            track=tracks[organelle],
            samples=samples,
            output_dir=output_dir,
            reference=reference,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            max_depth=max_depth,
            threads=threads,
            force=force,
        )
        results.append(result)
        all_command_rows.extend(command_rows)
    write_variant_calling_outputs(
        output_dir=output_dir,
        results=results,
        command_rows=all_command_rows,
        sample_count=len(samples),
        run_label=run_label,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        max_depth=max_depth,
    )
    return results


def labeled_output_name(name: str, run_label: str) -> str:
    if not run_label:
        return name
    return f"{run_label}.{name}"


def write_variant_calling_outputs(
    output_dir: Path,
    results: list[VariantCallResult],
    command_rows: list[dict[str, str]],
    sample_count: int,
    run_label: str,
    min_mapq: int,
    min_baseq: int,
    max_depth: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / labeled_output_name("commands.tsv", run_label),
        command_rows,
        ["organelle", "track_id", "step", "command"],
    )
    write_tsv(
        output_dir / labeled_output_name("variant_calling_summary.tsv", run_label),
        [
            {
                "organelle": result.organelle,
                "track_id": result.track_id,
                "sample_count": str(result.sample_count),
                "variant_records": str(result.variant_records),
                "raw_vcf_path": result.final_vcf_path.as_posix(),
                "raw_vcf_index_path": result.final_index_path.as_posix(),
                "log_path": result.log_path.as_posix(),
            }
            for result in results
        ],
        [
            "organelle",
            "track_id",
            "sample_count",
            "variant_records",
            "raw_vcf_path",
            "raw_vcf_index_path",
            "log_path",
        ],
    )
    write_report(
        output_dir / labeled_output_name("variant_calling_report.md", run_label),
        results=results,
        sample_count=sample_count,
        run_label=run_label,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        max_depth=max_depth,
    )


def write_report(
    path: Path,
    results: list[VariantCallResult],
    sample_count: int,
    run_label: str,
    min_mapq: int,
    min_baseq: int,
    max_depth: int,
) -> None:
    label = run_label or "full"
    lines = [
        "# Haploid Variant Calling",
        "",
        "This step calls raw haploid variants separately for cpDNA and mtDNA.",
        "Filtering and consensus generation happen in later steps.",
        "",
        "## Run",
        "",
        f"- Run label: `{label}`",
        f"- Samples called: {sample_count}",
        f"- Minimum mapping quality: {min_mapq}",
        f"- Minimum base quality: {min_baseq}",
        f"- Per-sample maximum pileup depth: {max_depth}",
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
                f"- Raw variant records: {result.variant_records}",
                f"- Raw VCF: `{result.final_vcf_path}`",
                f"- Index: `{result.final_index_path}`",
                f"- Log: `{result.log_path}`",
                "",
            ]
        )
    path.write_text("\n".join(lines))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Call raw haploid cpDNA/mtDNA variants."
    )
    parser.add_argument("--sample-table", type=Path, default=DEFAULT_SAMPLE_TABLE)
    parser.add_argument("--track-table", type=Path, default=DEFAULT_TRACK_TABLE)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--bam-dir", type=Path, default=DEFAULT_BAM_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--sample-limit", type=int)
    parser.add_argument("--sample-id", action="append", default=[])
    parser.add_argument("--run-label", default="")
    parser.add_argument("--min-mapq", type=int, default=DEFAULT_MIN_MAPQ)
    parser.add_argument("--min-baseq", type=int, default=DEFAULT_MIN_BASEQ)
    parser.add_argument("--max-depth", type=int, default=DEFAULT_MAX_DEPTH)
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS)
    parser.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    sample_ids = set(args.sample_id) if args.sample_id else None
    results = run_variant_calling(
        sample_table=args.sample_table,
        track_table=args.track_table,
        reference=args.reference,
        bam_dir=args.bam_dir,
        output_dir=args.output_dir,
        sample_limit=args.sample_limit,
        sample_ids=sample_ids,
        run_label=args.run_label,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        max_depth=args.max_depth,
        threads=args.threads,
        force=args.force,
    )
    for result in results:
        print(
            f"{result.organelle}: {result.variant_records} raw variant records "
            f"across {result.sample_count} samples"
        )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
