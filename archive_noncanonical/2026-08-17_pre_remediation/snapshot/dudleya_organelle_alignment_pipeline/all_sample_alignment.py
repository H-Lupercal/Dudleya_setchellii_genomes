"""Run all primary samples through cpDNA/mtDNA alignment and track-aware QC.

This stage scales the alignment/QC machinery to
all primary paired-end samples and applies the analysis tracks when
summarizing coverage. It does not call variants, build consensus FASTAs, or run
PCA/tree/Fst/admixture analyses.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

from dudleya_organelle_alignment_pipeline.pilot_alignment import (
    AlignmentError,
    AlignmentSample,
    build_organelle_summary_rows,
    build_sample_summary,
    count_fastq_records,
    fmt_float,
    outputs_are_ready,
    outputs_for_sample,
    parse_depth_file,
    parse_idxstats_file,
    read_alignment_samples,
    read_fai_lengths,
    require_reference_indexes,
    require_tools,
    run_alignment_commands,
    run_qc_commands,
    write_tsv,
)


DEFAULT_SAMPLE_TABLE = Path(
    "dudleya_organelle_alignment_pipeline/results/00_manifest/analysis_samples.tsv"
)
DEFAULT_TRACK_TABLE = Path(
    "dudleya_organelle_alignment_pipeline/results/05_analysis_masks/analysis_tracks.tsv"
)
DEFAULT_REFERENCE = Path(
    "dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment"
)


@dataclass(frozen=True)
class TrackRegion:
    track_id: str
    organelle: str
    purpose: str
    record: str
    start_1based: int
    end_1based: int
    name: str
    bed_path: Path
    step5_use: str
    notes: str

    @property
    def length_bp(self) -> int:
        return self.end_1based - self.start_1based + 1


@dataclass(frozen=True)
class TrackMetrics:
    track_id: str
    organelle: str
    purpose: str
    record: str
    region_count: int
    region_bp: int
    total_depth: int
    bases_ge_1x: int
    bases_ge_5x: int
    bases_ge_10x: int

    @property
    def mean_depth(self) -> float:
        if self.region_bp == 0:
            return 0.0
        return self.total_depth / self.region_bp

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
        if self.region_bp == 0:
            return 0.0
        return bases / self.region_bp


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_track_regions(track_table: Path) -> list[TrackRegion]:
    """Read the analysis-track manifest and its BED files."""

    regions: list[TrackRegion] = []
    for track in read_tsv(track_table):
        bed_path = Path(track["bed_path"])
        if not bed_path.exists():
            raise AlignmentError(
                f"Missing BED file for track {track['track_id']}: {bed_path}"
            )
        with bed_path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for fields in reader:
                if not fields or fields[0].startswith("#"):
                    continue
                if len(fields) < 3:
                    raise AlignmentError(f"Invalid BED row in {bed_path}: {fields}")
                record = fields[0]
                bed_start = int(fields[1])
                bed_end = int(fields[2])
                if bed_start < 0 or bed_end <= bed_start:
                    raise AlignmentError(f"Invalid BED interval in {bed_path}: {fields}")
                name = fields[3] if len(fields) > 3 else track["track_id"]
                regions.append(
                    TrackRegion(
                        track_id=track["track_id"],
                        organelle=track["organelle"],
                        purpose=track["purpose"],
                        record=record,
                        start_1based=bed_start + 1,
                        end_1based=bed_end,
                        name=name,
                        bed_path=bed_path,
                        step5_use=track.get("step5_use", ""),
                        notes=track.get("notes", ""),
                    )
                )

    if not regions:
        raise AlignmentError(f"No analysis-track regions found in {track_table}")
    return regions


def validate_track_regions(
    track_regions: list[TrackRegion],
    reference_lengths: dict[str, int],
) -> None:
    for region in track_regions:
        reference_length = reference_lengths.get(region.record)
        if reference_length is None:
            raise AlignmentError(
                f"Track {region.track_id} uses record {region.record}, "
                "which is absent from the reference FAI."
            )
        if region.end_1based > reference_length:
            raise AlignmentError(
                f"Track {region.track_id} interval {region.record}:"
                f"{region.start_1based}-{region.end_1based} exceeds "
                f"reference length {reference_length}."
            )


def parse_track_depth_file(
    depth_path: Path,
    track_regions: list[TrackRegion],
) -> dict[str, TrackMetrics]:
    """Summarize depth over the analysis-track intervals.

    Missing positions count as zero because each track's denominator is the
    total BED-defined region length, not only positions observed in the depth
    file.
    """

    counters = initialize_track_counters(track_regions)
    regions_by_record: dict[str, list[TrackRegion]] = {}
    for region in track_regions:
        regions_by_record.setdefault(region.record, []).append(region)

    with depth_path.open() as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            record = fields[0]
            position = int(fields[1])
            depth = int(fields[2])
            for region in regions_by_record.get(record, []):
                if region.start_1based <= position <= region.end_1based:
                    counter = counters[region.track_id]
                    counter["total_depth"] += depth
                    if depth >= 1:
                        counter["bases_ge_1x"] += 1
                    if depth >= 5:
                        counter["bases_ge_5x"] += 1
                    if depth >= 10:
                        counter["bases_ge_10x"] += 1

    return {
        track_id: TrackMetrics(
            track_id=track_id,
            organelle=values["organelle"],
            purpose=values["purpose"],
            record=values["record"],
            region_count=values["region_count"],
            region_bp=values["region_bp"],
            total_depth=values["total_depth"],
            bases_ge_1x=values["bases_ge_1x"],
            bases_ge_5x=values["bases_ge_5x"],
            bases_ge_10x=values["bases_ge_10x"],
        )
        for track_id, values in counters.items()
    }


def initialize_track_counters(track_regions: list[TrackRegion]) -> dict[str, dict]:
    counters: dict[str, dict] = {}
    for region in track_regions:
        if region.track_id not in counters:
            counters[region.track_id] = {
                "organelle": region.organelle,
                "purpose": region.purpose,
                "record": region.record,
                "region_count": 0,
                "region_bp": 0,
                "total_depth": 0,
                "bases_ge_1x": 0,
                "bases_ge_5x": 0,
                "bases_ge_10x": 0,
            }
        counters[region.track_id]["region_count"] += 1
        counters[region.track_id]["region_bp"] += region.length_bp
    return counters


def build_track_summary_rows(
    sample_id: str,
    row: dict[str, str],
    track_regions: list[TrackRegion],
    track_metrics: dict[str, TrackMetrics],
    depth_path: Path,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    track_beds = {
        region.track_id: region.bed_path.as_posix() for region in track_regions
    }
    for track_id in track_order(track_regions):
        metrics = track_metrics[track_id]
        rows.append(
            {
                "sample_id": sample_id,
                "batch": row.get("batch", ""),
                "species": row.get("species", ""),
                "popcode": row.get("popcode", ""),
                "track_id": metrics.track_id,
                "organelle": metrics.organelle,
                "purpose": metrics.purpose,
                "record": metrics.record,
                "region_count": str(metrics.region_count),
                "region_bp": str(metrics.region_bp),
                "total_depth": str(metrics.total_depth),
                "mean_depth": fmt_float(metrics.mean_depth),
                "bases_ge_1x": str(metrics.bases_ge_1x),
                "breadth_ge_1x": fmt_float(metrics.breadth_ge_1x),
                "bases_ge_5x": str(metrics.bases_ge_5x),
                "breadth_ge_5x": fmt_float(metrics.breadth_ge_5x),
                "bases_ge_10x": str(metrics.bases_ge_10x),
                "breadth_ge_10x": fmt_float(metrics.breadth_ge_10x),
                "bed_path": track_beds[track_id],
                "depth_path": depth_path.as_posix(),
            }
        )
    return rows


def track_order(track_regions: list[TrackRegion]) -> list[str]:
    ordered: list[str] = []
    for region in track_regions:
        if region.track_id not in ordered:
            ordered.append(region.track_id)
    return ordered


def run_all_sample_alignment(
    sample_table: Path,
    track_table: Path,
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
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    require_tools(("bwa", "samtools"))
    require_reference_indexes(reference_path)

    reference_lengths = read_fai_lengths(Path(f"{reference_path}.fai"))
    track_regions = read_track_regions(track_table)
    validate_track_regions(track_regions, reference_lengths)
    samples = read_alignment_samples(sample_table, sample_limit, sample_ids)
    if not samples:
        raise AlignmentError(f"No eligible primary samples found in {sample_table}")

    output_dir.mkdir(parents=True, exist_ok=True)
    sample_summaries: list[dict[str, str]] = []
    organelle_rows: list[dict[str, str]] = []
    track_rows: list[dict[str, str]] = []
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
                    "command": (
                        "outputs already present; pass --force or --refresh-qc "
                        "to regenerate"
                    ),
                }
            )

        input_read_records = 0
        if count_input_reads:
            input_read_records = count_fastq_records(sample.r1_path)
            input_read_records += count_fastq_records(sample.r2_path)

        mapped_counts = parse_idxstats_file(outputs.idxstats_path)
        depth_metrics = parse_depth_file(outputs.depth_path, reference_lengths)
        track_metrics = parse_track_depth_file(outputs.depth_path, track_regions)

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
        track_rows.extend(
            build_track_summary_rows(
                sample_id=sample.sample_id,
                row=sample.row,
                track_regions=track_regions,
                track_metrics=track_metrics,
                depth_path=outputs.depth_path,
            )
        )

        write_all_sample_outputs(
            output_dir=output_dir,
            sample_summaries=sample_summaries,
            organelle_rows=organelle_rows,
            track_rows=track_rows,
            command_rows=command_rows,
            reference_path=reference_path,
            sample_table=sample_table,
            track_table=track_table,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            completed_samples=len(sample_summaries),
            total_samples=len(samples),
        )

    return sample_summaries, organelle_rows, track_rows


def write_all_sample_outputs(
    output_dir: Path,
    sample_summaries: list[dict[str, str]],
    organelle_rows: list[dict[str, str]],
    track_rows: list[dict[str, str]],
    command_rows: list[dict[str, str]],
    reference_path: Path,
    sample_table: Path,
    track_table: Path,
    min_mapq: int,
    min_baseq: int,
    completed_samples: int,
    total_samples: int,
) -> None:
    if sample_summaries:
        write_tsv(
            output_dir / "all_sample_alignment_sample_summary.tsv",
            sample_summaries,
            list(sample_summaries[0].keys()),
        )
    if organelle_rows:
        write_tsv(
            output_dir / "all_sample_alignment_by_organelle.tsv",
            organelle_rows,
            list(organelle_rows[0].keys()),
        )
    if track_rows:
        write_tsv(
            output_dir / "all_sample_alignment_by_track.tsv",
            track_rows,
            list(track_rows[0].keys()),
        )
    if command_rows:
        write_tsv(
            output_dir / "commands.tsv",
            command_rows,
            ["sample_id", "step", "command"],
        )
    write_report(
        output_dir / "all_sample_alignment_report.md",
        sample_summaries=sample_summaries,
        organelle_rows=organelle_rows,
        track_rows=track_rows,
        reference_path=reference_path,
        sample_table=sample_table,
        track_table=track_table,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        completed_samples=completed_samples,
        total_samples=total_samples,
    )


def write_report(
    path: Path,
    sample_summaries: list[dict[str, str]],
    organelle_rows: list[dict[str, str]],
    track_rows: list[dict[str, str]],
    reference_path: Path,
    sample_table: Path,
    track_table: Path,
    min_mapq: int,
    min_baseq: int,
    completed_samples: int,
    total_samples: int,
) -> None:
    total_mapped = sum(
        int(row["total_organelle_mapped_reads"]) for row in sample_summaries
    )
    flagged = [
        row for row in sample_summaries if row["qc_notes"] != "pass_initial_mapping_screen"
    ]
    lines = [
        "# All-Sample Organelle Alignment",
        "",
        "This step maps every primary paired-end sample to the combined",
        "cpDNA/mtDNA reference and summarizes coverage using the",
        "analysis tracks. It does not call variants or build consensus FASTAs.",
        "",
        "## Inputs",
        "",
        f"- Sample table: `{sample_table}`",
        f"- Reference: `{reference_path}`",
        f"- Analysis tracks: `{track_table}`",
        f"- Minimum mapping quality retained in BAM/depth: `{min_mapq}`",
        f"- Minimum base quality used for depth: `{min_baseq}`",
        "",
        "## Progress",
        "",
        f"- Samples completed in this run: {completed_samples}",
        f"- Target samples for this invocation: {total_samples}",
        "",
        "## Summary",
        "",
        f"- Total cpDNA+mtDNA mapped read records: {total_mapped}",
        f"- Samples with initial QC notes: {len(flagged)}",
        "",
        "## Median Breadth At 1x By Organelle",
        "",
        f"- Chloroplast: {median_breadth(organelle_rows, 'chloroplast')}",
        f"- Mitochondria: {median_breadth(organelle_rows, 'mitochondria')}",
        "",
        "## Median Breadth At 1x By Analysis Track",
        "",
        *format_track_medians(track_rows),
        "",
        "## Outputs",
        "",
        "- `all_sample_alignment_sample_summary.tsv`: one row per sample.",
        "- `all_sample_alignment_by_organelle.tsv`: one row per sample and organelle.",
        "- `all_sample_alignment_by_track.tsv`: one row per sample and analysis track.",
        "- `commands.tsv`: external commands run plus reuse decisions.",
        "- `bam/`, `qc/`, and `logs/`: generated alignment artifacts ignored by git.",
        "",
        "Review this report before variant calling, consensus generation, PCA,",
        "tree building, Fst, or structure/admixture-style analyses.",
        "",
    ]
    if flagged:
        lines.extend(
            [
                "## Samples With QC Notes",
                "",
                *[
                    f"- `{row['sample_id']}`: {row['qc_notes']}"
                    for row in flagged[:50]
                ],
                "",
            ]
        )
    path.write_text("\n".join(lines))


def median_breadth(rows: list[dict[str, str]], key_value: str) -> str:
    values = sorted(
        float(row["breadth_ge_1x"])
        for row in rows
        if row.get("organelle") == key_value
    )
    if not values:
        return "NA"
    midpoint = len(values) // 2
    if len(values) % 2:
        return fmt_float(values[midpoint])
    return fmt_float((values[midpoint - 1] + values[midpoint]) / 2)


def format_track_medians(track_rows: list[dict[str, str]]) -> list[str]:
    if not track_rows:
        return ["- none"]
    by_track: dict[str, list[float]] = {}
    for row in track_rows:
        by_track.setdefault(row["track_id"], []).append(float(row["breadth_ge_1x"]))
    lines: list[str] = []
    for track_id in sorted(by_track):
        values = sorted(by_track[track_id])
        midpoint = len(values) // 2
        if len(values) % 2:
            median = values[midpoint]
        else:
            median = (values[midpoint - 1] + values[midpoint]) / 2
        lines.append(f"- {track_id}: {fmt_float(median)}")
    return lines


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run all-sample cpDNA/mtDNA alignment and track-aware QC."
    )
    parser.add_argument("--sample-table", type=Path, default=DEFAULT_SAMPLE_TABLE)
    parser.add_argument("--track-table", type=Path, default=DEFAULT_TRACK_TABLE)
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
        help="Skip FASTQ read counting; input mapping fractions will be reported as 0.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    sample_ids = set(args.sample_ids) if args.sample_ids else None
    sample_summaries, organelle_rows, track_rows = run_all_sample_alignment(
        sample_table=args.sample_table,
        track_table=args.track_table,
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
    print(f"All-sample alignment samples summarized: {len(sample_summaries)}")
    print(f"Organelle summary rows: {len(organelle_rows)}")
    print(f"Track summary rows: {len(track_rows)}")
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
