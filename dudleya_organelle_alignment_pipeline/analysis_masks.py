"""Define cpDNA and mtDNA analysis masks for all-sample processing.

This stage does not align reads or call variants.
It turns the cpDNA and mtDNA verification results into explicit BED tracks so
the all-sample run can apply the same scientific rules reproducibly.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

from dudleya_organelle_alignment_pipeline.prepare_reference_and_pilot import (
    EXPECTED_REFERENCE_LENGTHS,
)


DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/05_analysis_masks"
)
DEFAULT_CPDNA_REPEAT_INTERVALS = Path(
    "dudleya_organelle_alignment_pipeline/results/04_cpdna_investigation/"
    "cpdna_self_repeat_intervals.tsv"
)
DEFAULT_MTDNA_HIGH_MAPQ_INTERVALS = Path(
    "dudleya_organelle_alignment_pipeline/results/03_mtdna_investigation/"
    "mtdna_high_mapq_consensus_intervals.tsv"
)
DEFAULT_MTDNA_HIGH_CONFIDENCE_THRESHOLD = 12


class MaskDefinitionError(ValueError):
    """Raised when this stage source evidence cannot define a usable mask."""


@dataclass(frozen=True)
class Region:
    track_id: str
    organelle: str
    record: str
    start_1based: int
    end_1based: int
    name: str
    source: str
    note: str

    @property
    def length_bp(self) -> int:
        return self.end_1based - self.start_1based + 1


@dataclass(frozen=True)
class CpdnaTracks:
    full_coverage_regions: list[Region]
    ir_regions: list[Region]
    duplicate_ir_mask: list[Region]
    population_regions: list[Region]


@dataclass(frozen=True)
class MtdnaTracks:
    permissive_coverage_regions: list[Region]
    high_confidence_unique_regions: list[Region]


def interval_to_bed_fields(
    record: str,
    start_1based: int,
    end_1based: int,
    name: str,
) -> list[str]:
    """Convert 1-based inclusive coordinates to 0-based half-open BED fields."""

    validate_interval(start_1based, end_1based)
    return [record, str(start_1based - 1), str(end_1based), name]


def validate_interval(start_1based: int, end_1based: int) -> None:
    if start_1based < 1:
        raise MaskDefinitionError(f"Interval start must be >= 1: {start_1based}")
    if end_1based < start_1based:
        raise MaskDefinitionError(
            f"Interval end must be >= start: {start_1based}-{end_1based}"
        )


def validate_interval_within_reference(
    start_1based: int,
    end_1based: int,
    reference_length: int,
) -> None:
    validate_interval(start_1based, end_1based)
    if end_1based > reference_length:
        raise MaskDefinitionError(
            f"Interval {start_1based}-{end_1based} exceeds reference length "
            f"{reference_length}"
        )


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_cpdna_tracks(
    repeat_intervals_path: Path,
    reference_length: int = EXPECTED_REFERENCE_LENGTHS["chloroplast"],
) -> CpdnaTracks:
    """Build cpDNA tracks, keeping one IR copy for population analyses."""

    ir_intervals = read_major_cpdna_repeat_pair(repeat_intervals_path, reference_length)
    ir_regions = [
        Region(
            track_id="cpdna_ir_regions",
            organelle="cpDNA",
            record="chloroplast",
            start_1based=start,
            end_1based=end,
            name=f"cpdna_IR_copy_{index}",
            source=repeat_intervals_path.as_posix(),
            note="Major chloroplast inverted-repeat copy from cpDNA verification.",
        )
        for index, (start, end) in enumerate(ir_intervals, start=1)
    ]

    duplicate_start, duplicate_end = ir_intervals[-1]
    duplicate_mask = [
        Region(
            track_id="cpdna_duplicate_ir_mask",
            organelle="cpDNA",
            record="chloroplast",
            start_1based=duplicate_start,
            end_1based=duplicate_end,
            name="cpdna_duplicate_IR_copy_mask",
            source=repeat_intervals_path.as_posix(),
            note=(
                "Later IR copy excluded from the cpDNA population-genetic track "
                "so analyses keep only one IR copy."
            ),
        )
    ]

    return CpdnaTracks(
        full_coverage_regions=[
            Region(
                track_id="cpdna_full_coverage",
                organelle="cpDNA",
                record="chloroplast",
                start_1based=1,
                end_1based=reference_length,
                name="cpdna_full_reference",
                source="verified_chloroplast_reference",
                note="Full cpDNA reference for sample-level coverage QC.",
            )
        ],
        ir_regions=ir_regions,
        duplicate_ir_mask=duplicate_mask,
        population_regions=complement_regions(
            track_id="cpdna_population_sites",
            organelle="cpDNA",
            record="chloroplast",
            reference_length=reference_length,
            masked_regions=duplicate_mask,
            name_prefix="cpdna_population_single_IR_region",
            source=repeat_intervals_path.as_posix(),
            note=(
                "cpDNA region allowed for PCA, Fst, tree, and clustering inputs; "
                "one duplicate IR copy is excluded."
            ),
        ),
    )


def read_major_cpdna_repeat_pair(
    repeat_intervals_path: Path,
    reference_length: int,
) -> list[tuple[int, int]]:
    rows = read_tsv(repeat_intervals_path)
    if not rows:
        raise MaskDefinitionError(f"No cpDNA repeat rows found in {repeat_intervals_path}")

    max_length = max(int(row["length_bp"]) for row in rows)
    largest_rows = [row for row in rows if int(row["length_bp"]) == max_length]
    intervals: set[tuple[int, int]] = set()
    for row in largest_rows:
        for start_field, end_field in (
            ("query_start", "query_end"),
            ("match_start", "match_end"),
        ):
            start = int(row[start_field])
            end = int(row[end_field])
            validate_interval_within_reference(start, end, reference_length)
            intervals.add((start, end))

    ordered = sorted(intervals)
    if len(ordered) != 2:
        raise MaskDefinitionError(
            "Expected exactly two major cpDNA IR intervals from the largest "
            f"self-repeat evidence, found {len(ordered)}: {ordered}"
        )
    return ordered


def build_mtdna_tracks(
    high_mapq_intervals_path: Path,
    reference_length: int = EXPECTED_REFERENCE_LENGTHS["mitochondria"],
    high_confidence_threshold: int = DEFAULT_MTDNA_HIGH_CONFIDENCE_THRESHOLD,
) -> MtdnaTracks:
    """Build mtDNA tracks for permissive coverage QC and unique-site analyses."""

    high_confidence_regions = read_mtdna_high_confidence_regions(
        high_mapq_intervals_path,
        reference_length=reference_length,
        threshold=high_confidence_threshold,
    )
    return MtdnaTracks(
        permissive_coverage_regions=[
            Region(
                track_id="mtdna_permissive_coverage",
                organelle="mtDNA",
                record="mitochondria",
                start_1based=1,
                end_1based=reference_length,
                name="mtdna_permissive_coverage_full_reference",
                source="verified_mitochondrial_reference",
                note=(
                    "Full mtDNA reference used only for permissive MAPQ sample "
                    "presence, breadth, and depth QC."
                ),
            )
        ],
        high_confidence_unique_regions=high_confidence_regions,
    )


def read_mtdna_high_confidence_regions(
    high_mapq_intervals_path: Path,
    reference_length: int,
    threshold: int,
) -> list[Region]:
    rows = [
        row
        for row in read_tsv(high_mapq_intervals_path)
        if int(row["threshold_usable_samples"]) == threshold
    ]
    if not rows:
        raise MaskDefinitionError(
            "No mtDNA high-MAPQ consensus intervals found for threshold "
            f"{threshold} in {high_mapq_intervals_path}"
        )

    regions: list[Region] = []
    for index, row in enumerate(
        sorted(rows, key=lambda item: (int(item["rank"]), int(item["start"]))),
        start=1,
    ):
        start = int(row["start"])
        end = int(row["end"])
        validate_interval_within_reference(start, end, reference_length)
        regions.append(
            Region(
                track_id="mtdna_high_confidence_unique",
                organelle="mtDNA",
                record="mitochondria",
                start_1based=start,
                end_1based=end,
                name=f"mtdna_high_confidence_unique_region_{index}",
                source=high_mapq_intervals_path.as_posix(),
                note=(
                    "High-MAPQ consensus interval retained for mtDNA variant "
                    "calling and population-genetic analyses."
                ),
            )
        )
    return regions


def complement_regions(
    track_id: str,
    organelle: str,
    record: str,
    reference_length: int,
    masked_regions: list[Region],
    name_prefix: str,
    source: str,
    note: str,
) -> list[Region]:
    merged_masks = merge_coordinate_pairs(
        [(region.start_1based, region.end_1based) for region in masked_regions]
    )
    regions: list[Region] = []
    next_start = 1
    for mask_start, mask_end in merged_masks:
        if next_start < mask_start:
            regions.append(
                build_named_region(
                    track_id=track_id,
                    organelle=organelle,
                    record=record,
                    start_1based=next_start,
                    end_1based=mask_start - 1,
                    name_prefix=name_prefix,
                    index=len(regions) + 1,
                    source=source,
                    note=note,
                )
            )
        next_start = max(next_start, mask_end + 1)

    if next_start <= reference_length:
        regions.append(
            build_named_region(
                track_id=track_id,
                organelle=organelle,
                record=record,
                start_1based=next_start,
                end_1based=reference_length,
                name_prefix=name_prefix,
                index=len(regions) + 1,
                source=source,
                note=note,
            )
        )
    return regions


def build_named_region(
    track_id: str,
    organelle: str,
    record: str,
    start_1based: int,
    end_1based: int,
    name_prefix: str,
    index: int,
    source: str,
    note: str,
) -> Region:
    return Region(
        track_id=track_id,
        organelle=organelle,
        record=record,
        start_1based=start_1based,
        end_1based=end_1based,
        name=f"{name_prefix}_{index}",
        source=source,
        note=note,
    )


def merge_coordinate_pairs(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged: list[tuple[int, int]] = [ordered[0]]
    for start, end in ordered[1:]:
        previous_start, previous_end = merged[-1]
        if start <= previous_end + 1:
            merged[-1] = (previous_start, max(previous_end, end))
        else:
            merged.append((start, end))
    return merged


def generate_analysis_masks(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    cpdna_repeat_intervals_path: Path = DEFAULT_CPDNA_REPEAT_INTERVALS,
    mtdna_high_mapq_intervals_path: Path = DEFAULT_MTDNA_HIGH_MAPQ_INTERVALS,
    cpdna_length: int = EXPECTED_REFERENCE_LENGTHS["chloroplast"],
    mtdna_length: int = EXPECTED_REFERENCE_LENGTHS["mitochondria"],
    mtdna_high_confidence_threshold: int = DEFAULT_MTDNA_HIGH_CONFIDENCE_THRESHOLD,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    cpdna = build_cpdna_tracks(cpdna_repeat_intervals_path, cpdna_length)
    mtdna = build_mtdna_tracks(
        mtdna_high_mapq_intervals_path,
        reference_length=mtdna_length,
        high_confidence_threshold=mtdna_high_confidence_threshold,
    )

    track_files = {
        "cpdna_full_coverage": output_dir / "cpdna_full_coverage_regions.bed",
        "cpdna_ir_regions": output_dir / "cpdna_ir_regions.bed",
        "cpdna_duplicate_ir_mask": output_dir / "cpdna_duplicate_ir_copy_mask.bed",
        "cpdna_population_sites": output_dir / "cpdna_population_sites.bed",
        "mtdna_permissive_coverage": output_dir / "mtdna_permissive_coverage_regions.bed",
        "mtdna_high_confidence_unique": (
            output_dir / "mtdna_high_confidence_unique_regions.bed"
        ),
    }
    regions_by_track = {
        "cpdna_full_coverage": cpdna.full_coverage_regions,
        "cpdna_ir_regions": cpdna.ir_regions,
        "cpdna_duplicate_ir_mask": cpdna.duplicate_ir_mask,
        "cpdna_population_sites": cpdna.population_regions,
        "mtdna_permissive_coverage": mtdna.permissive_coverage_regions,
        "mtdna_high_confidence_unique": mtdna.high_confidence_unique_regions,
    }
    for track_id, regions in regions_by_track.items():
        write_bed(track_files[track_id], regions)

    write_region_manifest(output_dir / "analysis_regions.tsv", regions_by_track)
    write_track_manifest(
        output_dir / "analysis_tracks.tsv",
        track_files=track_files,
        cpdna_repeat_intervals_path=cpdna_repeat_intervals_path,
        mtdna_high_mapq_intervals_path=mtdna_high_mapq_intervals_path,
        mtdna_high_confidence_threshold=mtdna_high_confidence_threshold,
    )
    write_summary(
        output_dir / "mask_summary.md",
        cpdna=cpdna,
        mtdna=mtdna,
        cpdna_length=cpdna_length,
        mtdna_length=mtdna_length,
        cpdna_repeat_intervals_path=cpdna_repeat_intervals_path,
        mtdna_high_mapq_intervals_path=mtdna_high_mapq_intervals_path,
        mtdna_high_confidence_threshold=mtdna_high_confidence_threshold,
    )


def write_bed(path: Path, regions: list[Region]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for region in regions:
            writer.writerow(
                interval_to_bed_fields(
                    region.record,
                    region.start_1based,
                    region.end_1based,
                    region.name,
                )
            )


def write_region_manifest(
    path: Path,
    regions_by_track: dict[str, list[Region]],
) -> None:
    fieldnames = [
        "track_id",
        "organelle",
        "record",
        "name",
        "start_1based",
        "end_1based",
        "bed_start_0based",
        "bed_end_0based",
        "length_bp",
        "source",
        "note",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for track_id in regions_by_track:
            for region in regions_by_track[track_id]:
                writer.writerow(
                    {
                        "track_id": region.track_id,
                        "organelle": region.organelle,
                        "record": region.record,
                        "name": region.name,
                        "start_1based": region.start_1based,
                        "end_1based": region.end_1based,
                        "bed_start_0based": region.start_1based - 1,
                        "bed_end_0based": region.end_1based,
                        "length_bp": region.length_bp,
                        "source": region.source,
                        "note": region.note,
                    }
                )


def write_track_manifest(
    path: Path,
    track_files: dict[str, Path],
    cpdna_repeat_intervals_path: Path,
    mtdna_high_mapq_intervals_path: Path,
    mtdna_high_confidence_threshold: int,
) -> None:
    rows = [
        {
            "track_id": "cpdna_full_coverage",
            "organelle": "cpDNA",
            "purpose": "sample_qc_coverage",
            "bed_path": track_files["cpdna_full_coverage"].as_posix(),
            "coordinate_system": "BED 0-based half-open; source intervals 1-based inclusive",
            "source": "verified_chloroplast_reference",
            "step5_use": "Compute cpDNA permissive breadth/depth and mapped-read QC.",
            "notes": "Not the preferred site set for PCA, Fst, trees, or clustering.",
        },
        {
            "track_id": "cpdna_ir_regions",
            "organelle": "cpDNA",
            "purpose": "repeat_annotation",
            "bed_path": track_files["cpdna_ir_regions"].as_posix(),
            "coordinate_system": "BED 0-based half-open; source intervals 1-based inclusive",
            "source": cpdna_repeat_intervals_path.as_posix(),
            "step5_use": "Report cpDNA inverted-repeat coverage separately.",
            "notes": "Documents both major IR copies in the normalized chloroplast reference.",
        },
        {
            "track_id": "cpdna_duplicate_ir_mask",
            "organelle": "cpDNA",
            "purpose": "mask",
            "bed_path": track_files["cpdna_duplicate_ir_mask"].as_posix(),
            "coordinate_system": "BED 0-based half-open; source intervals 1-based inclusive",
            "source": cpdna_repeat_intervals_path.as_posix(),
            "step5_use": "Exclude this later IR copy from cpDNA population-site outputs.",
            "notes": "Implements the single-IR-copy strategy.",
        },
        {
            "track_id": "cpdna_population_sites",
            "organelle": "cpDNA",
            "purpose": "variant_calling_and_population_genetics",
            "bed_path": track_files["cpdna_population_sites"].as_posix(),
            "coordinate_system": "BED 0-based half-open; source intervals 1-based inclusive",
            "source": cpdna_repeat_intervals_path.as_posix(),
            "step5_use": "Allowed cpDNA sites for PCA, Fst, tree, and clustering inputs.",
            "notes": "Keeps one IR copy and excludes the duplicate IR copy.",
        },
        {
            "track_id": "mtdna_permissive_coverage",
            "organelle": "mtDNA",
            "purpose": "sample_qc_coverage",
            "bed_path": track_files["mtdna_permissive_coverage"].as_posix(),
            "coordinate_system": "BED 0-based half-open; source intervals 1-based inclusive",
            "source": "verified_mitochondrial_reference",
            "step5_use": "Compute mtDNA permissive MAPQ presence, breadth, and depth QC.",
            "notes": "Do not use this whole-reference track for mtDNA population-genetic variants.",
        },
        {
            "track_id": "mtdna_high_confidence_unique",
            "organelle": "mtDNA",
            "purpose": "variant_calling_and_population_genetics",
            "bed_path": track_files["mtdna_high_confidence_unique"].as_posix(),
            "coordinate_system": "BED 0-based half-open; source intervals 1-based inclusive",
            "source": mtdna_high_mapq_intervals_path.as_posix(),
            "step5_use": "Allowed mtDNA sites for variant calling, PCA, Fst, trees, and clustering.",
            "notes": (
                "Uses high-MAPQ consensus intervals supported by at least "
                f"{mtdna_high_confidence_threshold} usable pilot samples."
            ),
        },
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "track_id",
                "organelle",
                "purpose",
                "bed_path",
                "coordinate_system",
                "source",
                "step5_use",
                "notes",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    path: Path,
    cpdna: CpdnaTracks,
    mtdna: MtdnaTracks,
    cpdna_length: int,
    mtdna_length: int,
    cpdna_repeat_intervals_path: Path,
    mtdna_high_mapq_intervals_path: Path,
    mtdna_high_confidence_threshold: int,
) -> None:
    cpdna_population_bp = sum(region.length_bp for region in cpdna.population_regions)
    cpdna_masked_bp = sum(region.length_bp for region in cpdna.duplicate_ir_mask)
    mtdna_unique_bp = sum(
        region.length_bp for region in mtdna.high_confidence_unique_regions
    )
    lines = [
        "# Analysis Masks",
        "",
        "This step defines the cpDNA and mtDNA tracks that the all-sample run",
        "must use. It does not align reads, call variants, or create final",
        "population-genetic outputs.",
        "",
        "## Coordinate Systems",
        "",
        "- BED files are 0-based, half-open.",
        "- `analysis_regions.tsv` records the same intervals as 1-based inclusive",
        "  coordinates plus their BED coordinates.",
        "",
        "## cpDNA Tracks",
        "",
        f"- Source: `{cpdna_repeat_intervals_path}`",
        f"- Full cpDNA reference length: {cpdna_length} bp.",
        "- Strategy: keep one chloroplast IR copy for population-genetic outputs.",
        f"- Duplicate IR bases masked from cpDNA population sites: {cpdna_masked_bp}.",
        f"- cpDNA population-site bases retained: {cpdna_population_bp}.",
        "- Use `cpdna_full_coverage_regions.bed` for sample-level coverage QC.",
        "- Use `cpdna_population_sites.bed` for PCA, Fst, trees, and",
        "  admixture-style clustering inputs.",
        "",
        "## mtDNA Tracks",
        "",
        f"- Source: `{mtdna_high_mapq_intervals_path}`",
        f"- Full mtDNA reference length: {mtdna_length} bp.",
        "- Strategy: keep mtDNA in two tracks.",
        "- Use `mtdna_permissive_coverage_regions.bed` for sample-level",
        "  permissive MAPQ coverage QC.",
        "- Use `mtdna_high_confidence_unique_regions.bed` for mtDNA variant",
        "  calling and population genetics.",
        "- High-confidence mtDNA threshold: intervals supported by at least",
        f"  {mtdna_high_confidence_threshold} usable pilot samples.",
        f"- High-confidence mtDNA bases retained: {mtdna_unique_bp}.",
        "",
        "## Outputs",
        "",
        "- `analysis_tracks.tsv`: machine-readable track purpose and downstream use.",
        "- `analysis_regions.tsv`: machine-readable interval audit table.",
        "- `*.bed`: regions and masks consumed by later alignment/QC/variant steps.",
        "",
    ]
    path.write_text("\n".join(lines))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create cpDNA and mtDNA analysis masks."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Folder for the BED tracks and manifests.",
    )
    parser.add_argument(
        "--cpdna-repeat-intervals",
        type=Path,
        default=DEFAULT_CPDNA_REPEAT_INTERVALS,
        help="cpDNA self-repeat interval table from the cpDNA verification.",
    )
    parser.add_argument(
        "--mtdna-high-mapq-intervals",
        type=Path,
        default=DEFAULT_MTDNA_HIGH_MAPQ_INTERVALS,
        help="mtDNA high-MAPQ consensus interval table from mtDNA investigation.",
    )
    parser.add_argument(
        "--mtdna-high-confidence-threshold",
        type=int,
        default=DEFAULT_MTDNA_HIGH_CONFIDENCE_THRESHOLD,
        help="Usable-pilot-sample threshold for mtDNA high-confidence regions.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    generate_analysis_masks(
        output_dir=args.output_dir,
        cpdna_repeat_intervals_path=args.cpdna_repeat_intervals,
        mtdna_high_mapq_intervals_path=args.mtdna_high_mapq_intervals,
        mtdna_high_confidence_threshold=args.mtdna_high_confidence_threshold,
    )
    return 0
