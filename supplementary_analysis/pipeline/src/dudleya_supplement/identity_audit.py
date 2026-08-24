"""Executable conservative sample-identity audit."""

from __future__ import annotations

import hashlib
import shlex
import shutil
import statistics
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from .identity import MixedAlleleCall, classify_mixed_allele_samples, parse_structured_id
from .io import read_tsv, write_tsv


def _paths(value: str) -> list[str]:
    return [item for item in value.split(";") if item]


def _md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _revalidate_provider_row(root: Path, row: dict[str, str]) -> dict[str, str]:
    relative = row["resolved_source_path"]
    if not relative:
        return {
            **row,
            "supplementary_observed_md5": "",
            "supplementary_status": "DECLARED_MISSING_NOT_HASHABLE",
        }
    current = _md5(root / relative)
    return {
        **row,
        "supplementary_observed_md5": current,
        "supplementary_status": "PASS" if current == row["expected_md5"] else "FAIL",
    }


def _sketch_sample(root: Path, output: Path, sample: str, inputs: list[str]) -> None:
    if output.with_suffix(".msh").is_file():
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    files = " ".join(shlex.quote(value) for value in inputs)
    command = f"gzip -cd {files} | mash sketch -r -k 31 -s 100000 -I {shlex.quote(sample)} -o {shlex.quote(str(output))} -"
    subprocess.run(["bash", "-o", "pipefail", "-c", command], cwd=root, check=True)


def _sketch_samples(
    root: Path,
    work: Path,
    samples: list[dict[str, str]],
    *,
    workers: int = 4,
) -> list[Path]:
    if workers < 1 or workers > 4:
        raise ValueError("Identity sketch workers must be between 1 and 4")
    prefixes = [work / "samples" / row["sample_id"] for row in samples]
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="mash-sketch") as executor:
        futures = [
            executor.submit(
                _sketch_sample,
                root,
                prefix,
                row["sample_id"],
                _paths(row["r1_paths"]) + _paths(row["r2_paths"]),
            )
            for row, prefix in zip(samples, prefixes, strict=True)
        ]
        for future in futures:
            future.result()
    return [prefix.with_suffix(".msh") for prefix in prefixes]


def _mash_audit(root: Path, run_id: str, samples: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    if shutil.which("mash") is None:
        raise RuntimeError("Mash 2.3 is required for the approved raw-read identity audit")
    work = root / f"supplementary_analysis/work/{run_id}/identity/mash"
    complete = [row for row in samples if row["pair_status"] == "complete" and row["analysis_eligible"] == "yes"]
    sample_sketches = _sketch_samples(root, work, complete, workers=4)
    combined = work / "all_samples"
    if not combined.with_suffix(".msh").is_file():
        subprocess.run(["mash", "paste", str(combined), *map(str, sample_sketches)], cwd=root, check=True)
    distances = work / "all_pairs.dist.tsv"
    if not distances.is_file():
        with distances.open("w") as handle:
            subprocess.run(
                ["mash", "dist", str(combined.with_suffix(".msh")), str(combined.with_suffix(".msh"))], stdout=handle, check=True
            )
    controls = complete[: min(10, len(complete))]
    within_distances: list[float] = []
    for row in controls:
        r1, r2 = _paths(row["r1_paths"]), _paths(row["r2_paths"])
        left = work / "controls" / f"{row['sample_id']}.R1"
        right = work / "controls" / f"{row['sample_id']}.R2"
        _sketch_sample(root, left, f"{row['sample_id']}_R1", r1)
        _sketch_sample(root, right, f"{row['sample_id']}_R2", r2)
        result = subprocess.run(
            ["mash", "dist", str(left.with_suffix(".msh")), str(right.with_suffix(".msh"))], capture_output=True, text=True, check=True
        )
        within_distances.append(float(result.stdout.split("\t")[2]))
    different_controls = [
        (controls[index], controls[index + 1])
        for index in range(len(controls) - 1)
        if controls[index]["popcode"] != controls[index + 1]["popcode"]
    ]
    different_distances = []
    for left, right in different_controls:
        result = subprocess.run(
            [
                "mash",
                "dist",
                str((work / "samples" / left["sample_id"]).with_suffix(".msh")),
                str((work / "samples" / right["sample_id"]).with_suffix(".msh")),
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        different_distances.append(float(result.stdout.split("\t")[2]))
    within_cutoff = max(within_distances)
    different_floor = min(different_distances) if different_distances else -1.0
    cutoff = min(within_cutoff, different_floor) if different_floor >= 0 else within_cutoff
    rows = []
    seen: set[tuple[str, str]] = set()
    for line in distances.read_text().splitlines():
        left, right, distance, p_value, shared = line.split("\t")
        left, right = Path(left).name, Path(right).name
        if left == right:
            continue
        pair = tuple(sorted((left, right)))
        if pair in seen:
            continue
        seen.add(pair)
        value = float(distance)
        rows.append(
            {
                "sample_1": pair[0],
                "sample_2": pair[1],
                "mash_distance": f"{value:.12g}",
                "p_value": p_value,
                "shared_hashes": shared,
                "calibrated_candidate": "yes" if value <= cutoff else "no",
                "interpretation": "suspected_only" if value <= cutoff else "not_flagged_negative_screen_not_proof",
            }
        )
    calibration = [
        {
            "control_type": "within_library_R1_vs_R2",
            "count": len(within_distances),
            "minimum": min(within_distances),
            "median": statistics.median(within_distances),
            "maximum": max(within_distances),
            "candidate_cutoff": cutoff,
        },
        {
            "control_type": "different_library",
            "count": len(different_distances),
            "minimum": min(different_distances) if different_distances else "NA",
            "median": statistics.median(different_distances) if different_distances else "NA",
            "maximum": max(different_distances) if different_distances else "NA",
            "candidate_cutoff": cutoff,
        },
    ]
    return rows, calibration


def _mixed_allele_audit(root: Path) -> list[dict[str, object]]:
    calls: dict[str, list[MixedAlleleCall]] = defaultdict(list)
    for organelle in ("chloroplast", "mitochondria"):
        vcf = root / f"canonical_publication/results/variants/publication-20260817/{organelle}.high_confidence_variant_sites.vcf.gz"
        samples = subprocess.run(["bcftools", "query", "-l", str(vcf)], capture_output=True, text=True, check=True).stdout.splitlines()
        process = subprocess.Popen(["bcftools", "query", "-f", "%POS[\t%DP:%AD]\n", str(vcf)], stdout=subprocess.PIPE, text=True)
        assert process.stdout is not None
        for line in process.stdout:
            fields = line.rstrip().split("\t")[1:]
            for sample, value in zip(samples, fields, strict=True):
                depth_text, alleles = value.split(":", 1)
                if depth_text == "." or alleles == ".":
                    continue
                depths = [int(item) for item in alleles.split(",") if item != "."]
                if len(depths) >= 2:
                    calls[sample].append(MixedAlleleCall(int(depth_text), depths[0], sum(depths[1:])))
        if process.wait() != 0:
            raise RuntimeError(f"bcftools mixed-allele query failed for {organelle}")
    return [row.__dict__ for row in classify_mixed_allele_samples(calls)]


def run_identity_audit(root: Path, run_id: str) -> list[Path]:
    samples = read_tsv(root / "canonical_publication/metadata/samples/samples.tsv")
    provider = read_tsv(root / "canonical_publication/provenance/manifests/publication-20260817.provider_md5_validation.tsv")
    output = root / f"supplementary_analysis/results/verification/{run_id}/identity"
    structured = []
    for row in samples:
        for read_path in _paths(row["r1_paths"]) + _paths(row["r2_paths"]):
            structured.append({"sample_id": row["sample_id"], "read_path": read_path, **parse_structured_id(Path(read_path).name)})
    path_to_samples: dict[str, set[str]] = defaultdict(set)
    for row in samples:
        for value in _paths(row["r1_paths"]) + _paths(row["r2_paths"]):
            path_to_samples[value].add(row["sample_id"])
    revalidated = []
    exact_groups: dict[str, list[str]] = defaultdict(list)
    for row in provider:
        result = _revalidate_provider_row(root, row)
        revalidated.append(result)
        if row["resolved_source_path"]:
            exact_groups[result["supplementary_observed_md5"]].append(row["resolved_source_path"])
    exact = [
        {
            "observed_md5": digest,
            "file_count": len(paths),
            "sample_count": len(set().union(*(path_to_samples[path] for path in paths))),
            "sample_ids": ";".join(sorted(set().union(*(path_to_samples[path] for path in paths)))),
            "paths": ";".join(sorted(paths)),
            "status": "confirmed_cross_sample_exact_duplicate"
            if len(set().union(*(path_to_samples[path] for path in paths))) > 1
            else "same_sample_or_non_sample_duplicate",
        }
        for digest, paths in exact_groups.items()
        if len(set(paths)) > 1
    ]
    mash_rows, calibration = _mash_audit(root, run_id, samples)
    mixed = _mixed_allele_audit(root)
    index = [
        {
            "status": "untestable",
            "reason": "index sequences, sample sheet, demultiplexing metrics, and unexpected index combinations are unavailable",
            "S_number_interpretation": "demultiplexed sample number; not an index sequence",
        }
    ]
    paths = {
        "provider": output / "provider_md5_revalidation.tsv",
        "structured": output / "structured_identifier_audit.tsv",
        "exact": output / "exact_duplicate_files.tsv",
        "mash": output / "read_sketch_pairwise.tsv",
        "calibration": output / "read_sketch_calibration.tsv",
        "mixed": output / "mixed_allele_screen.tsv",
        "index": output / "index_hopping_assessability.tsv",
    }
    write_tsv(paths["provider"], revalidated, list(revalidated[0]), root)
    write_tsv(paths["structured"], structured, ["sample_id", "read_path", "plate_well", "specimen", "demultiplex_sample", "lane"], root)
    write_tsv(paths["exact"], exact, ["observed_md5", "file_count", "sample_count", "sample_ids", "paths", "status"], root)
    write_tsv(paths["mash"], mash_rows, list(mash_rows[0]), root)
    write_tsv(paths["calibration"], calibration, list(calibration[0]), root)
    write_tsv(paths["mixed"], mixed, list(mixed[0]), root)
    write_tsv(paths["index"], index, list(index[0]), root)
    outcomes = []
    suspected = {row["sample_1"] for row in mash_rows if row["calibrated_candidate"] == "yes"} | {
        row["sample_2"] for row in mash_rows if row["calibrated_candidate"] == "yes"
    }
    mixed_suspected = {str(row["sample_id"]) for row in mixed if row["status"] == "suspected"}
    confirmed = {
        sample for row in exact if row["status"] == "confirmed_cross_sample_exact_duplicate" for sample in str(row["sample_ids"]).split(";")
    }
    for row in samples:
        sample = row["sample_id"]
        status = (
            "confirmed_duplicate"
            if sample in confirmed
            else "suspected"
            if sample in suspected or sample in mixed_suspected
            else "unresolved"
        )
        outcomes.append(
            {
                "sample_id": sample,
                "outcome": status,
                "action": "correct_or_exclude_and_invalidate"
                if status.startswith("confirmed")
                else "with_and_without_sensitivity"
                if status == "suspected"
                else "retain_with_limitation_no_replication_or_clonality_claim",
            }
        )
    outcome_path = output / "sample_identity_outcomes.tsv"
    write_tsv(outcome_path, outcomes, ["sample_id", "outcome", "action"], root)
    return [*paths.values(), outcome_path]
