#!/usr/bin/env python3
"""Compute organelle-specific breadth, eligibility, and the mtDNA site mask."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import subprocess
import tomllib
from array import array
from dataclasses import asdict
from pathlib import Path

from organelle_pipeline.paths import repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)
from organelle_pipeline.qc import (
    SampleBreadth,
    build_depth_command,
    select_organelle_samples,
    summarize_depths,
    summarize_fastp_report,
    summarize_masked_depths,
    support_intervals,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def reference_lengths(fai: Path) -> dict[str, int]:
    result = {}
    for line in fai.read_text().splitlines():
        fields = line.split("\t")
        result[fields[0]] = int(fields[1])
    return result


def bed_inclusion_mask(path: Path, record: str, length: int) -> bytearray:
    mask = bytearray(length)
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        observed_record, start, end, *_ = line.split("\t")
        if observed_record == record:
            for index in range(int(start), int(end)):
                mask[index] = 1
    return mask


def read_depths(
    bam: Path,
    record: str,
    length: int,
    minimum_mapping_quality: int,
    minimum_base_quality: int,
) -> array:
    depths = array("I", [0]) * length
    process = subprocess.Popen(
        build_depth_command(
            bam,
            record,
            minimum_mapping_quality=minimum_mapping_quality,
            minimum_base_quality=minimum_base_quality,
        ),
        stdout=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    for line in process.stdout:
        _, position, depth = line.rstrip().split("\t")
        depths[int(position) - 1] = int(depth)
    if process.wait() != 0:
        raise RuntimeError(f"samtools depth failed for {bam} {record}")
    return depths


def write_rows(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    sample_manifest = root / config["paths"]["sample_manifest"]
    reference = root / "canonical_publication/references/selected/organelle_combined.fa"
    lengths = reference_lengths(Path(f"{reference}.fai"))
    complete_samples = [row for row in read_tsv(sample_manifest) if row["analysis_eligible"] == "yes" and row["pair_status"] == "complete"]
    bam_dir = root / "canonical_publication/work" / args.run_id / "mapping"
    mapping_state_dir = root / "canonical_publication/provenance/runs" / args.run_id / "mapping"
    fastp_dir = root / "canonical_publication/results/qc" / args.run_id / "fastp"
    output_dir = root / "canonical_publication/results/qc" / args.run_id
    metadata_dir = root / "canonical_publication/metadata/qc" / args.run_id
    mask_dir = root / "canonical_publication/references/masks" / args.run_id
    state_path = root / "canonical_publication/provenance/runs" / args.run_id / "qc.json"
    for directory in (output_dir, metadata_dir, mask_dir, state_path.parent):
        directory.mkdir(parents=True, exist_ok=True)
    reference_state_path = root / "canonical_publication/provenance/runs" / args.run_id / "references.json"
    reference_state = json.loads(reference_state_path.read_text())
    mapping_completion_path = root / "canonical_publication/provenance/runs" / args.run_id / "mapping_complete.json"
    if not mapping_completion_path.is_file():
        raise RuntimeError("QC requires finalized dependency-complete mapping provenance")
    mapping_completion = json.loads(mapping_completion_path.read_text())
    if mapping_completion.get("status") != "complete":
        raise RuntimeError("Mapping provenance completion state is not successful")
    validate_saved_outputs(root, mapping_completion)
    cp_eligibility_mask = root / "canonical_publication/references/masks/chloroplast_unique_eligibility_sites.bed"
    mt_repeat_mask = root / "canonical_publication/references/masks/mitochondria_repeat_mask.bed"
    mt_annotations = root / "canonical_publication/references/annotations/mitochondria.projected.tsv"
    upstream = {
        "references": reference_state["fingerprint"]["digest"],
        "mapping_completion": mapping_completion["fingerprint"]["digest"],
    }
    declared_inputs = {
        **runtime_provenance(root, {"samtools": ("samtools", "--version")}),
        sample_manifest.relative_to(root).as_posix(): sha256_file(sample_manifest),
        reference.relative_to(root).as_posix(): sha256_file(reference),
        cp_eligibility_mask.relative_to(root).as_posix(): sha256_file(cp_eligibility_mask),
        mt_repeat_mask.relative_to(root).as_posix(): sha256_file(mt_repeat_mask),
        mt_annotations.relative_to(root).as_posix(): sha256_file(mt_annotations),
        config_path.relative_to(root).as_posix(): sha256_file(config_path),
        mapping_completion_path.relative_to(root).as_posix(): sha256_file(mapping_completion_path),
    }
    for row in complete_samples:
        sample_id = row["sample_id"]
        bam = bam_dir / f"{sample_id}.organelle.bam"
        fastp_json = fastp_dir / f"{sample_id}.fastp.json"
        state = mapping_state_dir / f"{sample_id}.json"
        if not bam.exists() or not fastp_json.exists() or not state.exists():
            raise RuntimeError(f"Missing canonical mapping or provenance for {sample_id}")
        mapping_state = json.loads(state.read_text())
        bam_key = bam.relative_to(root).as_posix()
        bam_hash = sha256_file(bam)
        if mapping_state["outputs"].get(bam_key) != bam_hash:
            raise RuntimeError(f"Mapping output checksum mismatch for {sample_id}")
        fastp_key = fastp_json.relative_to(root).as_posix()
        fastp_hash = sha256_file(fastp_json)
        if mapping_state["outputs"].get(fastp_key) != fastp_hash:
            raise RuntimeError(f"fastp output checksum mismatch for {sample_id}")
        declared_inputs[bam_key] = bam_hash
        declared_inputs[fastp_key] = fastp_hash
        upstream[sample_id] = mapping_state["fingerprint"]["digest"]
    minimum_mapping_quality = int(config["mapping"]["minimum_mapping_quality"])
    minimum_base_quality = int(config["mapping"]["minimum_base_quality"])
    commands = [
        " ".join(
            build_depth_command(
                "SAMPLE.bam",
                record,
                minimum_mapping_quality=minimum_mapping_quality,
                minimum_base_quality=minimum_base_quality,
            )
        )
        for record in ("chloroplast", "mitochondria")
    ]
    fingerprint = build_stage_fingerprint_from_hashes("qc", declared_inputs, upstream, commands)
    if args.resume and state_path.exists():
        saved = json.loads(state_path.read_text())
        validate_resume(saved["fingerprint"]["digest"], fingerprint)
        for path, digest in saved["outputs"].items():
            if sha256_file(root / path) != digest:
                raise RuntimeError(f"QC output checksum mismatch: {path}")
        print("resume-valid qc")
        return 0
    if state_path.exists():
        raise RuntimeError("QC state already exists; use --resume or a new run ID")
    expected_qc_outputs = [
        output_dir / "sample_breadth.tsv",
        output_dir / "read_preprocessing_summary.tsv",
        metadata_dir / "organelle_eligibility.tsv",
        metadata_dir / "chloroplast_samples.tsv",
        metadata_dir / "mitochondria_samples.tsv",
        metadata_dir / "shared_samples.tsv",
        mask_dir / "mitochondria_high_confidence_sites.bed",
        root / "canonical_publication/references/evidence" / args.run_id / "read_backed_reference_validation.tsv",
        root / "canonical_publication/references/evidence" / args.run_id / "mitochondria_high_confidence_annotation_overlap.tsv",
        root / "canonical_publication/references/evidence" / args.run_id / "mitochondria_repeat_read_support.tsv",
    ]
    if any(path.exists() for path in expected_qc_outputs):
        raise RuntimeError("Existing unvalidated QC output; preserve it and use a new run ID")

    breadth_thresholds = tuple(int(value) for value in config["qc"]["breadth_depths"])
    eligibility_depth = int(config["qc"]["eligibility_depth"])
    if config["qc"]["breadth_denominator"] != "organelle_unique_mappability_mask":
        raise RuntimeError("Unsupported QC breadth denominator")
    if eligibility_depth not in breadth_thresholds:
        raise RuntimeError("QC eligibility depth must be one of the reported breadth depths")
    cp_analysis_mask = bed_inclusion_mask(
        cp_eligibility_mask,
        "chloroplast",
        lengths["chloroplast"],
    )
    mt_unique_mask = bytearray([1]) * lengths["mitochondria"]
    for line in mt_repeat_mask.read_text().splitlines():
        if not line.strip():
            continue
        _, start, end, *_ = line.split("\t")
        for index in range(int(start), int(end)):
            mt_unique_mask[index] = 0
    mt_repeat_sites = bytearray(not bool(value) for value in mt_unique_mask)

    # Retain one byte per mitochondrial position and sample until eligibility is
    # known.  High-confidence support must be calculated from the independently
    # eligible mitochondrial set, not boosted or diluted by ineligible samples.
    mt_callable_by_sample: dict[str, bytearray] = {}
    mt_repeat_depth_by_sample: dict[str, tuple[float, float]] = {}
    cp_window = min(1000, lengths["chloroplast"] // 4)
    mt_window = min(1000, lengths["mitochondria"] // 4)
    boundary_depths: dict[str, dict[str, tuple[float, float, float]]] = {
        "chloroplast": {},
        "mitochondria": {},
    }
    summaries = []
    breadth_rows = []
    preprocessing_rows = []
    for row in complete_samples:
        sample_id = row["sample_id"]
        bam = bam_dir / f"{sample_id}.organelle.bam"
        preprocessing = summarize_fastp_report(json.loads((fastp_dir / f"{sample_id}.fastp.json").read_text()))
        preprocessing_rows.append(
            {
                "sample_id": sample_id,
                **{key: f"{value:.8f}" if isinstance(value, float) else value for key, value in asdict(preprocessing).items()},
            }
        )
        cp_depths = read_depths(
            bam,
            "chloroplast",
            lengths["chloroplast"],
            minimum_mapping_quality,
            minimum_base_quality,
        )
        mt_depths = read_depths(
            bam,
            "mitochondria",
            lengths["mitochondria"],
            minimum_mapping_quality,
            minimum_base_quality,
        )
        cp_full = summarize_depths(cp_depths, breadth_thresholds)
        mt_full = summarize_depths(mt_depths, breadth_thresholds)
        cp = summarize_masked_depths(cp_depths, cp_analysis_mask, breadth_thresholds)
        mt = summarize_masked_depths(mt_depths, mt_unique_mask, breadth_thresholds)
        mt_repeat = summarize_masked_depths(mt_depths, mt_repeat_sites, breadth_thresholds)
        mt_repeat_depth_by_sample[sample_id] = (mt_repeat.mean_depth, mt.mean_depth)
        boundary_depths["chloroplast"][sample_id] = (
            sum(cp_depths[:cp_window]) / cp_window,
            sum(cp_depths[-cp_window:]) / cp_window,
            sum(cp_depths[cp_window:-cp_window]) / (lengths["chloroplast"] - 2 * cp_window),
        )
        boundary_depths["mitochondria"][sample_id] = (
            sum(mt_depths[:mt_window]) / mt_window,
            sum(mt_depths[-mt_window:]) / mt_window,
            sum(mt_depths[mt_window:-mt_window]) / (lengths["mitochondria"] - 2 * mt_window),
        )
        mt_callable_by_sample[sample_id] = bytearray(
            depth >= eligibility_depth and bool(mt_unique_mask[index]) for index, depth in enumerate(mt_depths)
        )
        summaries.append(
            SampleBreadth(
                sample_id,
                cp_dp5=cp.breadth[eligibility_depth],
                mt_dp5=mt.breadth[eligibility_depth],
            )
        )
        breadth_rows.append(
            {
                "sample_id": sample_id,
                "cp_full_reference_mean_depth": f"{cp_full.mean_depth:.6f}",
                **{f"cp_full_reference_breadth_dp{depth}": f"{cp_full.breadth[depth]:.8f}" for depth in breadth_thresholds},
                "cp_unique_sites_mean_depth": f"{cp.mean_depth:.6f}",
                **{f"cp_unique_sites_breadth_dp{depth}": f"{cp.breadth[depth]:.8f}" for depth in breadth_thresholds},
                "mt_full_reference_mean_depth": f"{mt_full.mean_depth:.6f}",
                **{f"mt_full_reference_breadth_dp{depth}": f"{mt_full.breadth[depth]:.8f}" for depth in breadth_thresholds},
                "mt_unique_sites_mean_depth": f"{mt.mean_depth:.6f}",
                **{f"mt_unique_sites_breadth_dp{depth}": f"{mt.breadth[depth]:.8f}" for depth in breadth_thresholds},
                "mt_repeat_sites_mean_depth": f"{mt_repeat.mean_depth:.6f}",
                **{f"mt_repeat_sites_breadth_dp{depth}": f"{mt_repeat.breadth[depth]:.8f}" for depth in breadth_thresholds},
            }
        )
        print(f"qc {sample_id}", flush=True)
    minimum_breadth = float(config["qc"]["minimum_breadth"])
    selected = select_organelle_samples(summaries, minimum_breadth)
    if not selected.cp or not selected.mt:
        raise RuntimeError("No samples passed one or both organelle-specific QC eligibility rules")
    cp_set, mt_set, shared_set = set(selected.cp), set(selected.mt), set(selected.shared)
    mt_support = array("H", [0]) * lengths["mitochondria"]
    for sample_id in selected.mt:
        for index, callable_position in enumerate(mt_callable_by_sample[sample_id]):
            mt_support[index] += callable_position
    eligibility_rows = []
    for row in complete_samples:
        sample_id = row["sample_id"]
        eligibility_rows.append(
            {
                "sample_id": sample_id,
                "include_chloroplast": "yes" if sample_id in cp_set else "no",
                "include_mitochondria": "yes" if sample_id in mt_set else "no",
                "include_shared": "yes" if sample_id in shared_set else "no",
                "rule": f"organelle_unique_mappability_breadth_dp{eligibility_depth}>={minimum_breadth:.2f}",
                "manual_override": "no",
            }
        )
    breadth_path = output_dir / "sample_breadth.tsv"
    preprocessing_path = output_dir / "read_preprocessing_summary.tsv"
    eligibility_path = metadata_dir / "organelle_eligibility.tsv"
    write_rows(breadth_path, breadth_rows, list(breadth_rows[0]))
    write_rows(
        preprocessing_path,
        preprocessing_rows,
        list(preprocessing_rows[0]),
    )
    write_rows(
        eligibility_path,
        eligibility_rows,
        [
            "sample_id",
            "include_chloroplast",
            "include_mitochondria",
            "include_shared",
            "rule",
            "manual_override",
        ],
    )
    metadata = {row["sample_id"]: row for row in complete_samples}
    for label, sample_ids in (
        ("chloroplast", selected.cp),
        ("mitochondria", selected.mt),
        ("shared", selected.shared),
    ):
        rows = [
            {
                "sample_id": sample_id,
                "popcode": metadata[sample_id]["popcode"],
                "species": metadata[sample_id]["species"],
                "population_name": metadata[sample_id]["population_name"],
            }
            for sample_id in sample_ids
        ]
        write_rows(
            metadata_dir / f"{label}_samples.tsv",
            rows,
            ["sample_id", "popcode", "species", "population_name"],
        )

    intervals = support_intervals(
        mt_support,
        sample_count=len(selected.mt),
        minimum_fraction=float(config["qc"]["mitochondria_high_confidence_minimum_sample_fraction"]),
        minimum_length=int(config["qc"]["mitochondria_high_confidence_minimum_interval_length"]),
    )
    if not intervals:
        raise RuntimeError("No mitochondrial interval passed the read-backed high-confidence rule")
    mt_mask = mask_dir / "mitochondria_high_confidence_sites.bed"
    with mt_mask.open("w") as handle:
        for number, (start, end) in enumerate(intervals, 1):
            handle.write(f"mitochondria\t{start}\t{end}\thigh_confidence_{number}\n")
    evidence_dir = root / "canonical_publication/references/evidence" / args.run_id
    evidence_dir.mkdir(parents=True, exist_ok=True)
    read_validation = evidence_dir / "read_backed_reference_validation.tsv"
    annotation_overlap = evidence_dir / "mitochondria_high_confidence_annotation_overlap.tsv"
    repeat_support = evidence_dir / "mitochondria_repeat_read_support.tsv"

    repeat_ratios = []
    with repeat_support.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample_id",
                "repeat_mask_mean_depth",
                "unique_sites_mean_depth",
                "repeat_to_unique_mean_depth_ratio",
                "interpretation",
            ]
        )
        for sample_id in selected.mt:
            repeat_mean, unique_mean = mt_repeat_depth_by_sample[sample_id]
            ratio = repeat_mean / unique_mean if unique_mean else float("nan")
            if ratio == ratio:
                repeat_ratios.append(ratio)
            writer.writerow(
                [
                    sample_id,
                    f"{repeat_mean:.8f}",
                    f"{unique_mean:.8f}",
                    f"{ratio:.8f}",
                    "MAPQ/base-quality-filtered depth diagnostic; repeat coordinates remain excluded from inference",
                ]
            )

    annotation_rows = []
    for feature in read_tsv(mt_annotations):
        feature_start = int(feature["start_1based"]) - 1
        feature_end = int(feature["end_1based"])
        feature_length = feature_end - feature_start
        overlap_bases = sum(
            max(0, min(feature_end, interval_end) - max(feature_start, interval_start)) for interval_start, interval_end in intervals
        )
        annotation_rows.append(
            {
                "feature_id": feature["feature_id"],
                "feature_type": feature["feature_type"],
                "gene": feature["gene"],
                "start_1based": feature["start_1based"],
                "end_1based": feature["end_1based"],
                "feature_length": feature_length,
                "high_confidence_overlap_bases": overlap_bases,
                "high_confidence_overlap_fraction": f"{overlap_bases / feature_length:.8f}",
                "overlap_status": ("full" if overlap_bases == feature_length else ("partial" if overlap_bases else "none")),
                "annotation_status": feature["status"],
            }
        )
    write_rows(
        annotation_overlap,
        annotation_rows,
        [
            "feature_id",
            "feature_type",
            "gene",
            "start_1based",
            "end_1based",
            "feature_length",
            "high_confidence_overlap_bases",
            "high_confidence_overlap_fraction",
            "overlap_status",
            "annotation_status",
        ],
    )

    def eligible_boundary_summary(
        organelle: str,
        sample_ids: tuple[str, ...],
    ) -> tuple[float, float, float, float]:
        values = [boundary_depths[organelle][sample_id] for sample_id in sample_ids]
        start_mean = statistics.fmean(value[0] for value in values)
        end_mean = statistics.fmean(value[1] for value in values)
        interior_mean = statistics.fmean(value[2] for value in values)
        ratios = [min(start, end) / interior for start, end, interior in values if interior > 0]
        return start_mean, end_mean, interior_mean, statistics.median(ratios) if ratios else float("nan")

    cp_start_mean, cp_end_mean, cp_interior_mean, cp_boundary_ratio = eligible_boundary_summary("chloroplast", selected.cp)
    mt_start_mean, mt_end_mean, mt_interior_mean, mt_boundary_ratio = eligible_boundary_summary("mitochondria", selected.mt)
    mt_high_confidence_bases = sum(end - start for start, end in intervals)
    with read_validation.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value", "interpretation"])
        writer.writerow(
            [
                "chloroplast_start_1kb_eligible_sample_mean_depth",
                f"{cp_start_mean:.8f}",
                "circular boundary flank; equal weight per chloroplast-eligible sample",
            ]
        )
        writer.writerow(
            [
                "chloroplast_unique_eligibility_bases",
                sum(cp_analysis_mask),
                "both inverted-repeat copies excluded from breadth denominator",
            ]
        )
        writer.writerow(
            [
                "chloroplast_end_1kb_eligible_sample_mean_depth",
                f"{cp_end_mean:.8f}",
                "circular boundary flank; equal weight per chloroplast-eligible sample",
            ]
        )
        writer.writerow(
            [
                "chloroplast_interior_eligible_sample_mean_depth",
                f"{cp_interior_mean:.8f}",
                "reference interior; equal weight per chloroplast-eligible sample",
            ]
        )
        writer.writerow(
            [
                "chloroplast_median_sample_boundary_to_interior_ratio",
                f"{cp_boundary_ratio:.8f}",
                "median among chloroplast-eligible samples; low values flag unsupported circularization boundaries",
            ]
        )
        writer.writerow(
            [
                "mitochondria_start_1kb_eligible_sample_mean_depth",
                f"{mt_start_mean:.8f}",
                "assembly boundary flank; equal weight per mitochondrial-eligible sample",
            ]
        )
        writer.writerow(
            [
                "mitochondria_unique_eligibility_bases",
                sum(mt_unique_mask),
                "both copies of qualifying self-repeats excluded from breadth denominator",
            ]
        )
        writer.writerow(
            [
                "mitochondria_end_1kb_eligible_sample_mean_depth",
                f"{mt_end_mean:.8f}",
                "assembly boundary flank; equal weight per mitochondrial-eligible sample",
            ]
        )
        writer.writerow(
            [
                "mitochondria_interior_eligible_sample_mean_depth",
                f"{mt_interior_mean:.8f}",
                "reference interior; equal weight per mitochondrial-eligible sample",
            ]
        )
        writer.writerow(
            [
                "mitochondria_median_sample_boundary_to_interior_ratio",
                f"{mt_boundary_ratio:.8f}",
                "median among mitochondrial-eligible samples; interpret with repeat mask",
            ]
        )
        writer.writerow(
            [
                "mitochondria_high_confidence_unique_bases",
                mt_high_confidence_bases,
                f"DP>={eligibility_depth} in "
                f">={float(config['qc']['mitochondria_high_confidence_minimum_sample_fraction']):.0%} "
                "of mitochondrial-eligible samples after self-repeat exclusion",
            ]
        )
        writer.writerow(
            [
                "mitochondria_median_sample_repeat_to_unique_mean_depth_ratio",
                f"{statistics.median(repeat_ratios):.8f}" if repeat_ratios else "nan",
                "eligible-sample read-backed repeat diagnostic; repeat coordinates remain excluded regardless of depth",
            ]
        )
    outputs = {
        path.relative_to(root).as_posix(): sha256_file(path)
        for path in [
            breadth_path,
            preprocessing_path,
            eligibility_path,
            metadata_dir / "chloroplast_samples.tsv",
            metadata_dir / "mitochondria_samples.tsv",
            metadata_dir / "shared_samples.tsv",
            mt_mask,
            read_validation,
            annotation_overlap,
            repeat_support,
        ]
    }
    state_path.write_text(
        json.dumps(
            {
                "status": "complete",
                "fingerprint": asdict(fingerprint),
                "sample_counts": {
                    "chloroplast": len(selected.cp),
                    "mitochondria": len(selected.mt),
                    "shared": len(selected.shared),
                },
                "outputs": outputs,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(
        f"eligible cp={len(selected.cp)} mt={len(selected.mt)} shared={len(selected.shared)}; mt high-confidence intervals={len(intervals)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
