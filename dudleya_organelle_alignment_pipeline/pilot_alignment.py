"""Run pilot cpDNA/mtDNA read mapping and summarize organelle signal.

This is step 3 of the pipeline. It aligns only the representative pilot samples
chosen in step 2, keeps mapped reads against the combined cpDNA/mtDNA reference,
and writes mapping/depth summaries for review before any all-sample run.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_REFERENCE = Path(
    "dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa"
)
DEFAULT_PILOT_TABLE = Path(
    "dudleya_organelle_alignment_pipeline/results/01_reference_pilot/pilot_samples.tsv"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/02_pilot_alignment"
)
ORGANELLES = ("chloroplast", "mitochondria")
LOW_MAPPED_READS_THRESHOLD = 100
LOW_BREADTH_THRESHOLD = 0.50


class AlignmentError(RuntimeError):
    """Raised when the pilot alignment step cannot safely continue."""


@dataclass(frozen=True)
class AlignmentSample:
    sample_id: str
    row: dict[str, str]
    r1_path: Path
    r2_path: Path


@dataclass(frozen=True)
class AlignmentOutputs:
    bam_path: Path
    bam_index_path: Path
    flagstat_path: Path
    idxstats_path: Path
    depth_path: Path
    log_path: Path


@dataclass(frozen=True)
class OrganelleMetrics:
    organelle: str
    reference_length: int
    total_depth: int
    bases_ge_1x: int
    bases_ge_5x: int
    bases_ge_10x: int

    @property
    def mean_depth(self) -> float:
        if self.reference_length == 0:
            return 0.0
        return self.total_depth / self.reference_length

    @property
    def breadth_ge_1x(self) -> float:
        return self._breadth(self.bases_ge_1x)

    @property
    def breadth_ge_5x(self) -> float:
        return self._breadth(self.bases_ge_5x)

    @property
    def breadth_ge_10x(self) -> float:
        return self._breadth(self.bases_ge_10x)

    def _breadth(self, bases: int) -> float:
        if self.reference_length == 0:
            return 0.0
        return bases / self.reference_length


def safe_sample_name(sample_id: str) -> str:
    """Make a sample ID safe for filenames while preserving useful IDs."""

    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", sample_id.strip())
    cleaned = cleaned.strip("_")
    return cleaned or "sample"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_path_field(value: str) -> list[Path]:
    return [Path(item) for item in value.split(";") if item]


def read_alignment_samples(
    table_path: Path,
    sample_limit: int | None = None,
    sample_ids: set[str] | None = None,
) -> list[AlignmentSample]:
    """Read pilot rows and keep only complete paired-end primary-analysis samples."""

    samples: list[AlignmentSample] = []
    for row in read_tsv(table_path):
        sample_id = row.get("sample_id", "")
        if sample_ids is not None and sample_id not in sample_ids:
            continue
        if row.get("analysis_status") != "include_primary_paired_end":
            continue
        if row.get("pair_status") != "complete":
            continue

        r1_paths = split_path_field(row.get("r1_paths", ""))
        r2_paths = split_path_field(row.get("r2_paths", ""))
        if not r1_paths or not r2_paths:
            continue
        if len(r1_paths) != 1 or len(r2_paths) != 1:
            raise AlignmentError(
                "Step 3 currently expects exactly one R1 and one R2 FASTQ per "
                f"pilot sample. Review sample {sample_id} before alignment."
            )
        samples.append(
            AlignmentSample(
                sample_id=sample_id,
                row=row,
                r1_path=r1_paths[0],
                r2_path=r2_paths[0],
            )
        )
        if sample_limit is not None and len(samples) >= sample_limit:
            break
    return samples


def count_fastq_records(path: Path) -> int:
    """Count FASTQ records by streaming newline counts in text or gzip files."""

    opener = gzip.open if path.suffix == ".gz" else open
    line_count = 0
    with opener(path, "rb") as handle:  # type: ignore[arg-type]
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            line_count += chunk.count(b"\n")
    if line_count % 4 != 0:
        raise AlignmentError(f"FASTQ line count is not divisible by 4: {path}")
    return line_count // 4


def read_fai_lengths(fai_path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with fai_path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def parse_idxstats_file(path: Path) -> dict[str, int]:
    mapped_counts: dict[str, int] = {}
    with path.open() as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] == "*":
                continue
            mapped_counts[fields[0]] = int(fields[2])
    return mapped_counts


def parse_depth_file(path: Path, reference_lengths: dict[str, int]) -> dict[str, OrganelleMetrics]:
    counters = {
        organelle: {
            "total_depth": 0,
            "bases_ge_1x": 0,
            "bases_ge_5x": 0,
            "bases_ge_10x": 0,
        }
        for organelle in reference_lengths
    }
    with path.open() as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            organelle = fields[0]
            if organelle not in counters:
                continue
            depth = int(fields[2])
            counters[organelle]["total_depth"] += depth
            if depth >= 1:
                counters[organelle]["bases_ge_1x"] += 1
            if depth >= 5:
                counters[organelle]["bases_ge_5x"] += 1
            if depth >= 10:
                counters[organelle]["bases_ge_10x"] += 1

    return {
        organelle: OrganelleMetrics(
            organelle=organelle,
            reference_length=reference_lengths[organelle],
            total_depth=values["total_depth"],
            bases_ge_1x=values["bases_ge_1x"],
            bases_ge_5x=values["bases_ge_5x"],
            bases_ge_10x=values["bases_ge_10x"],
        )
        for organelle, values in counters.items()
    }


def fmt_float(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}f}"


def build_sample_summary(
    sample_id: str,
    row: dict[str, str],
    mapped_counts: dict[str, int],
    depth_metrics: dict[str, OrganelleMetrics],
    input_read_records: int | None = None,
) -> dict[str, str]:
    cp_mapped = mapped_counts.get("chloroplast", 0)
    mt_mapped = mapped_counts.get("mitochondria", 0)
    total_mapped = cp_mapped + mt_mapped
    cp_fraction = cp_mapped / total_mapped if total_mapped else 0.0
    mt_fraction = mt_mapped / total_mapped if total_mapped else 0.0
    input_mapping_fraction = (
        total_mapped / input_read_records if input_read_records else 0.0
    )
    cp_metrics = depth_metrics.get(
        "chloroplast", OrganelleMetrics("chloroplast", 0, 0, 0, 0, 0)
    )
    mt_metrics = depth_metrics.get(
        "mitochondria", OrganelleMetrics("mitochondria", 0, 0, 0, 0, 0)
    )

    notes: list[str] = []
    if total_mapped == 0:
        notes.append("no_organelle_mapped_reads")
    if cp_mapped < LOW_MAPPED_READS_THRESHOLD:
        notes.append("low_chloroplast_mapped_reads")
    if mt_mapped < LOW_MAPPED_READS_THRESHOLD:
        notes.append("low_mitochondria_mapped_reads")
    if cp_metrics.breadth_ge_1x < LOW_BREADTH_THRESHOLD:
        notes.append("low_chloroplast_breadth_ge_1x")
    if mt_metrics.breadth_ge_1x < LOW_BREADTH_THRESHOLD:
        notes.append("low_mitochondria_breadth_ge_1x")

    return {
        "sample_id": sample_id,
        "batch": row.get("batch", ""),
        "naming_profile": row.get("naming_profile", ""),
        "species": row.get("species", ""),
        "popcode": row.get("popcode", ""),
        "input_read_records": str(input_read_records or 0),
        "total_organelle_mapped_reads": str(total_mapped),
        "input_organelle_mapping_fraction": fmt_float(input_mapping_fraction),
        "chloroplast_mapped_reads": str(cp_mapped),
        "mitochondria_mapped_reads": str(mt_mapped),
        "chloroplast_fraction_of_organelle_mapped": fmt_float(cp_fraction),
        "mitochondria_fraction_of_organelle_mapped": fmt_float(mt_fraction),
        "chloroplast_mean_depth": fmt_float(cp_metrics.mean_depth),
        "mitochondria_mean_depth": fmt_float(mt_metrics.mean_depth),
        "chloroplast_breadth_ge_1x": fmt_float(cp_metrics.breadth_ge_1x),
        "mitochondria_breadth_ge_1x": fmt_float(mt_metrics.breadth_ge_1x),
        "chloroplast_breadth_ge_5x": fmt_float(cp_metrics.breadth_ge_5x),
        "mitochondria_breadth_ge_5x": fmt_float(mt_metrics.breadth_ge_5x),
        "chloroplast_breadth_ge_10x": fmt_float(cp_metrics.breadth_ge_10x),
        "mitochondria_breadth_ge_10x": fmt_float(mt_metrics.breadth_ge_10x),
        "qc_notes": ";".join(notes) if notes else "pass_initial_mapping_screen",
    }


def build_organelle_summary_rows(
    sample: AlignmentSample,
    mapped_counts: dict[str, int],
    depth_metrics: dict[str, OrganelleMetrics],
    outputs: AlignmentOutputs,
    input_read_records: int,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for organelle in ORGANELLES:
        metrics = depth_metrics[organelle]
        rows.append(
            {
                "sample_id": sample.sample_id,
                "batch": sample.row.get("batch", ""),
                "species": sample.row.get("species", ""),
                "popcode": sample.row.get("popcode", ""),
                "organelle": organelle,
                "reference_length": str(metrics.reference_length),
                "input_read_records": str(input_read_records),
                "mapped_reads": str(mapped_counts.get(organelle, 0)),
                "mean_depth": fmt_float(metrics.mean_depth),
                "total_depth": str(metrics.total_depth),
                "bases_ge_1x": str(metrics.bases_ge_1x),
                "breadth_ge_1x": fmt_float(metrics.breadth_ge_1x),
                "bases_ge_5x": str(metrics.bases_ge_5x),
                "breadth_ge_5x": fmt_float(metrics.breadth_ge_5x),
                "bases_ge_10x": str(metrics.bases_ge_10x),
                "breadth_ge_10x": fmt_float(metrics.breadth_ge_10x),
                "bam_path": outputs.bam_path.as_posix(),
                "idxstats_path": outputs.idxstats_path.as_posix(),
                "depth_path": outputs.depth_path.as_posix(),
            }
        )
    return rows


def outputs_for_sample(output_dir: Path, sample_id: str) -> AlignmentOutputs:
    stem = safe_sample_name(sample_id)
    bam_path = output_dir / "bam" / f"{stem}.organelle.sorted.bam"
    return AlignmentOutputs(
        bam_path=bam_path,
        bam_index_path=Path(f"{bam_path}.bai"),
        flagstat_path=output_dir / "qc" / f"{stem}.flagstat.txt",
        idxstats_path=output_dir / "qc" / f"{stem}.idxstats.tsv",
        depth_path=output_dir / "qc" / f"{stem}.depth.tsv",
        log_path=output_dir / "logs" / f"{stem}.alignment.log",
    )


def outputs_are_ready(outputs: AlignmentOutputs) -> bool:
    return all(
        path.exists()
        for path in (
            outputs.bam_path,
            outputs.bam_index_path,
            outputs.flagstat_path,
            outputs.idxstats_path,
            outputs.depth_path,
        )
    )


def require_tools(tools: Iterable[str]) -> None:
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        raise AlignmentError(
            "Missing required alignment tools on PATH: " + ", ".join(sorted(missing))
        )


def require_reference_indexes(reference_path: Path) -> None:
    required = [
        Path(f"{reference_path}.fai"),
        Path(f"{reference_path}.amb"),
        Path(f"{reference_path}.ann"),
        Path(f"{reference_path}.bwt"),
        Path(f"{reference_path}.pac"),
        Path(f"{reference_path}.sa"),
    ]
    missing = [path.as_posix() for path in required if not path.exists()]
    if missing:
        raise AlignmentError(
            "Missing reference index files. Re-run Step 2 before alignment: "
            + "; ".join(missing)
        )


def build_depth_command(
    bam_path: Path,
    min_mapq: int,
    min_baseq: int,
) -> list[str]:
    """Build samtools depth command with explicit samtools flag semantics.

    In `samtools depth`, `-q` means minimum base quality and `-Q` means
    minimum mapping quality.
    """

    return [
        "samtools",
        "depth",
        "-aa",
        "-q",
        str(min_baseq),
        "-Q",
        str(min_mapq),
        bam_path.as_posix(),
    ]


def shlex_join(args: list[str | Path]) -> str:
    return shlex.join(str(arg) for arg in args)


def run_alignment_commands(
    sample: AlignmentSample,
    reference_path: Path,
    outputs: AlignmentOutputs,
    threads: int,
    min_mapq: int,
    min_baseq: int,
) -> list[dict[str, str]]:
    """Run bwa/samtools for one sample and return command-audit rows."""

    for path in (
        outputs.bam_path.parent,
        outputs.flagstat_path.parent,
        outputs.log_path.parent,
    ):
        path.mkdir(parents=True, exist_ok=True)

    tmp_bam = outputs.bam_path.with_suffix(".tmp.bam")
    align_command = [
        "bwa",
        "mem",
        "-t",
        str(threads),
        reference_path.as_posix(),
        sample.r1_path.as_posix(),
        sample.r2_path.as_posix(),
    ]
    view_command = [
        "samtools",
        "view",
        "-@",
        str(max(1, threads)),
        "-b",
        "-F",
        "4",
        "-q",
        str(min_mapq),
        "-",
    ]
    sort_command = [
        "samtools",
        "sort",
        "-@",
        str(max(1, threads)),
        "-o",
        tmp_bam.as_posix(),
        "-",
    ]
    command_rows = [
        {
            "sample_id": sample.sample_id,
            "step": "align_filter_sort",
            "command": " | ".join(
                shlex_join(command)
                for command in (align_command, view_command, sort_command)
            ),
        }
    ]

    with outputs.log_path.open("w") as log_handle:
        bwa_proc = subprocess.Popen(
            align_command,
            stdout=subprocess.PIPE,
            stderr=log_handle,
        )
        assert bwa_proc.stdout is not None
        view_proc = subprocess.Popen(
            view_command,
            stdin=bwa_proc.stdout,
            stdout=subprocess.PIPE,
            stderr=log_handle,
        )
        bwa_proc.stdout.close()
        assert view_proc.stdout is not None
        sort_proc = subprocess.Popen(
            sort_command,
            stdin=view_proc.stdout,
            stdout=subprocess.DEVNULL,
            stderr=log_handle,
        )
        view_proc.stdout.close()

        sort_return = sort_proc.wait()
        view_return = view_proc.wait()
        bwa_return = bwa_proc.wait()

    if bwa_return or view_return or sort_return:
        raise AlignmentError(
            f"Alignment failed for {sample.sample_id}; see {outputs.log_path}"
        )

    tmp_bam.replace(outputs.bam_path)
    command_rows.extend(
        run_qc_commands(
            sample=sample,
            outputs=outputs,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
        )
    )

    return command_rows


def run_qc_commands(
    sample: AlignmentSample,
    outputs: AlignmentOutputs,
    min_mapq: int,
    min_baseq: int,
) -> list[dict[str, str]]:
    """Create BAM index and mapping/depth QC files from an existing BAM."""

    outputs.flagstat_path.parent.mkdir(parents=True, exist_ok=True)
    outputs.log_path.parent.mkdir(parents=True, exist_ok=True)
    post_commands: list[tuple[str, list[str], Path | None]] = [
        (
            "index",
            ["samtools", "index", outputs.bam_path.as_posix()],
            None,
        ),
        (
            "flagstat",
            ["samtools", "flagstat", outputs.bam_path.as_posix()],
            outputs.flagstat_path,
        ),
        (
            "idxstats",
            ["samtools", "idxstats", outputs.bam_path.as_posix()],
            outputs.idxstats_path,
        ),
        (
            "depth",
            build_depth_command(outputs.bam_path, min_mapq, min_baseq),
            outputs.depth_path,
        ),
    ]
    command_rows: list[dict[str, str]] = []
    for step, command, stdout_path in post_commands:
        run_logged_command(command, stdout_path, outputs.log_path)
        command_rows.append(
            {
                "sample_id": sample.sample_id,
                "step": step,
                "command": shlex_join(command),
            }
        )

    return command_rows


def run_logged_command(command: list[str], stdout_path: Path | None, log_path: Path) -> None:
    with log_path.open("a") as log_handle:
        if stdout_path is None:
            result = subprocess.run(command, stderr=log_handle, check=False)
        else:
            with stdout_path.open("w") as stdout_handle:
                result = subprocess.run(
                    command,
                    stdout=stdout_handle,
                    stderr=log_handle,
                    check=False,
                )
    if result.returncode != 0:
        raise AlignmentError(f"Command failed: {shlex_join(command)}; see {log_path}")


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    sample_summaries: list[dict[str, str]],
    organelle_rows: list[dict[str, str]],
    reference_path: Path,
    min_mapq: int,
    min_baseq: int,
) -> None:
    sample_count = len(sample_summaries)
    total_mapped = sum(
        int(row["total_organelle_mapped_reads"]) for row in sample_summaries
    )
    flagged = [
        row for row in sample_summaries if row["qc_notes"] != "pass_initial_mapping_screen"
    ]
    cp_rows = [row for row in organelle_rows if row["organelle"] == "chloroplast"]
    mt_rows = [row for row in organelle_rows if row["organelle"] == "mitochondria"]
    lines = [
        "# Step 3 Pilot Organelle Alignment",
        "",
        "This step aligns the representative pilot samples to the combined",
        "cpDNA/mtDNA reference and summarizes organelle mapping signal. It does",
        "not call variants, make consensus FASTAs, or build final alignments.",
        "",
        "## Inputs",
        "",
        f"- Reference: `{reference_path}`",
        f"- Minimum mapping quality retained in BAM/depth: `{min_mapq}`",
        f"- Minimum base quality used for depth: `{min_baseq}`",
        "",
        "## Summary",
        "",
        f"- Samples attempted: {sample_count}",
        f"- Total cpDNA+mtDNA mapped read records: {total_mapped}",
        f"- Samples with initial QC notes: {len(flagged)}",
        "",
        "## Median Breadth At 1x",
        "",
        f"- Chloroplast: {median_breadth(cp_rows)}",
        f"- Mitochondria: {median_breadth(mt_rows)}",
        "",
        "## Outputs",
        "",
        "- `pilot_alignment_sample_summary.tsv`: one row per sample.",
        "- `pilot_alignment_by_organelle.tsv`: one row per sample and organelle.",
        "- `commands.tsv`: external commands run plus any reuse decisions.",
        "- `bam/`: filtered, sorted, indexed organelle BAM files.",
        "- `qc/`: per-sample flagstat, idxstats, and depth files.",
        "",
        "Review the sample and organelle summaries before scaling to all primary",
        "paired-end samples.",
        "",
    ]
    if flagged:
        lines.extend(
            [
                "## Samples With QC Notes",
                "",
                *[
                    f"- `{row['sample_id']}`: {row['qc_notes']}"
                    for row in flagged[:25]
                ],
                "",
            ]
        )
    path.write_text("\n".join(lines))


def median_breadth(rows: list[dict[str, str]]) -> str:
    values = sorted(float(row["breadth_ge_1x"]) for row in rows)
    if not values:
        return "NA"
    midpoint = len(values) // 2
    if len(values) % 2:
        return fmt_float(values[midpoint])
    return fmt_float((values[midpoint - 1] + values[midpoint]) / 2)


def run_pilot_alignment(
    pilot_table: Path,
    reference_path: Path,
    output_dir: Path,
    threads: int,
    min_mapq: int,
    min_baseq: int,
    sample_limit: int | None = None,
    sample_ids: set[str] | None = None,
    force: bool = False,
    refresh_qc: bool = False,
    count_input_reads: bool = True,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    require_tools(("bwa", "samtools"))
    require_reference_indexes(reference_path)

    reference_lengths = read_fai_lengths(Path(f"{reference_path}.fai"))
    samples = read_alignment_samples(pilot_table, sample_limit, sample_ids)
    if not samples:
        raise AlignmentError(f"No eligible pilot samples found in {pilot_table}")

    sample_summaries: list[dict[str, str]] = []
    organelle_rows: list[dict[str, str]] = []
    command_rows: list[dict[str, str]] = []

    for index, sample in enumerate(samples, start=1):
        print(f"[{index}/{len(samples)}] {sample.sample_id}", flush=True)
        outputs = outputs_for_sample(output_dir, sample.sample_id)
        if force or not outputs.bam_path.exists():
            command_rows.extend(
                run_alignment_commands(
                    sample=sample,
                    reference_path=reference_path,
                    outputs=outputs,
                    threads=threads,
                    min_mapq=min_mapq,
                    min_baseq=min_baseq,
                )
            )
        elif refresh_qc or not outputs_are_ready(outputs):
            command_rows.extend(
                run_qc_commands(
                    sample=sample,
                    outputs=outputs,
                    min_mapq=min_mapq,
                    min_baseq=min_baseq,
                )
            )
        else:
            command_rows.append(
                {
                    "sample_id": sample.sample_id,
                    "step": "reuse_existing_outputs",
                    "command": "outputs already present; pass --force or --refresh-qc to regenerate",
                }
            )

        input_read_records = 0
        if count_input_reads:
            input_read_records = count_fastq_records(sample.r1_path)
            input_read_records += count_fastq_records(sample.r2_path)

        mapped_counts = parse_idxstats_file(outputs.idxstats_path)
        depth_metrics = parse_depth_file(outputs.depth_path, reference_lengths)
        sample_summaries.append(
            build_sample_summary(
                sample_id=sample.sample_id,
                row=sample.row,
                mapped_counts=mapped_counts,
                depth_metrics=depth_metrics,
                input_read_records=input_read_records,
            )
        )
        organelle_rows.extend(
            build_organelle_summary_rows(
                sample=sample,
                mapped_counts=mapped_counts,
                depth_metrics=depth_metrics,
                outputs=outputs,
                input_read_records=input_read_records,
            )
        )

    write_tsv(
        output_dir / "pilot_alignment_sample_summary.tsv",
        sample_summaries,
        list(sample_summaries[0].keys()),
    )
    write_tsv(
        output_dir / "pilot_alignment_by_organelle.tsv",
        organelle_rows,
        list(organelle_rows[0].keys()),
    )
    write_tsv(
        output_dir / "commands.tsv",
        command_rows,
        ["sample_id", "step", "command"],
    )
    write_report(
        output_dir / "pilot_alignment_report.md",
        sample_summaries,
        organelle_rows,
        reference_path,
        min_mapq,
        min_baseq,
    )
    return sample_summaries, organelle_rows


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run Step 3 pilot cpDNA/mtDNA alignment and mapping QC."
    )
    parser.add_argument("--pilot-table", type=Path, default=DEFAULT_PILOT_TABLE)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--sample-limit", type=int)
    parser.add_argument(
        "--sample-id",
        action="append",
        dest="sample_ids",
        help="Restrict to one sample ID. May be provided multiple times.",
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--refresh-qc",
        action="store_true",
        help="Reuse existing BAMs but regenerate BAM indexes, flagstat, idxstats, and depth files.",
    )
    parser.add_argument(
        "--skip-input-read-counts",
        action="store_true",
        help="Skip FASTQ read counting; mapping fractions will be reported as 0.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    sample_ids = set(args.sample_ids) if args.sample_ids else None
    sample_summaries, organelle_rows = run_pilot_alignment(
        pilot_table=args.pilot_table,
        reference_path=args.reference,
        output_dir=args.output_dir,
        threads=args.threads,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        sample_limit=args.sample_limit,
        sample_ids=sample_ids,
        force=args.force,
        refresh_qc=args.refresh_qc,
        count_input_reads=not args.skip_input_read_counts,
    )
    print(f"Pilot samples summarized: {len(sample_summaries)}")
    print(f"Organelle summary rows: {len(organelle_rows)}")
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
