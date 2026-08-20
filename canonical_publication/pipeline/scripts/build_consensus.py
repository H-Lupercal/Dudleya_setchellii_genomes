#!/usr/bin/env python3
"""Build per-organelle callable alignments with uncertain sites encoded as N."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tomllib
from dataclasses import asdict
from pathlib import Path

from organelle_pipeline.consensus import analysis_mask_length, build_callable_sequence, reference_concordance
from organelle_pipeline.paths import repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
)
from organelle_pipeline.references import read_single_fasta, write_fasta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def read_samples(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        return [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")]


def read_mask_intervals(path: Path, record: str) -> list[tuple[int, int]]:
    intervals = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3 or fields[0] != record:
            raise ValueError(f"Unexpected analysis-mask row in {path}: {line}")
        intervals.append((int(fields[1]), int(fields[2])))
    return intervals


def accepted_variant_positions(vcf: Path) -> set[int]:
    completed = subprocess.run(
        ["bcftools", "query", "-f", "%POS\n", str(vcf)],
        capture_output=True,
        text=True,
        check=True,
    )
    return {int(value) for value in completed.stdout.splitlines() if value}


def sample_rows(bcf: Path, sample: str) -> list[tuple[int, str, str, str]]:
    process = subprocess.Popen(
        [
            "bcftools",
            "query",
            "-s",
            sample,
            "-f",
            "%POS\t%REF\t%ALT[\t%GT]\n",
            str(bcf),
        ],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    rows = []
    for line in process.stdout:
        position, ref, alt, genotype = line.rstrip().split("\t")
        rows.append((int(position), ref, alt, genotype))
    if process.wait() != 0:
        raise RuntimeError(f"bcftools query failed for {sample}")
    return rows


def sample_invariant_rows(
    bcf: Path,
    sample: str,
) -> list[tuple[int, str, str, int, tuple[int, ...]]]:
    """Read raw per-position depth and likelihoods before bcftools call."""

    process = subprocess.Popen(
        [
            "bcftools",
            "query",
            "-s",
            sample,
            "-f",
            "%POS\t%REF\t%ALT[\t%DP\t%PL]\n",
            str(bcf),
        ],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    rows = []
    for line in process.stdout:
        position, ref, alt, depth_text, likelihood_text = line.rstrip().split("\t")
        depth = int(depth_text) if depth_text not in {"", "."} else -1
        likelihoods = tuple(int(value) for value in likelihood_text.split(",")) if likelihood_text not in {"", "."} else ()
        rows.append((int(position), ref, alt, depth, likelihoods))
    if process.wait() != 0:
        raise RuntimeError(f"bcftools likelihood query failed for {sample}")
    return rows


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    minimum_depth = int(config["variants"]["minimum_depth"])
    minimum_genotype_quality = int(config["variants"]["minimum_genotype_quality"])
    references = dict(
        read_single_fasta(
            root / "canonical_publication/references/selected/organelle_combined.fa",
            expected_records=2,
        )
    )
    result_dir = root / "canonical_publication/results/alignments" / args.run_id
    variant_dir = root / "canonical_publication/results/variants" / args.run_id
    work_dir = root / "canonical_publication/work" / args.run_id / "variants"
    metadata_dir = root / "canonical_publication/metadata/qc" / args.run_id
    state_dir = root / "canonical_publication/provenance/runs" / args.run_id / "consensus"
    result_dir.mkdir(parents=True, exist_ok=True)
    state_dir.mkdir(parents=True, exist_ok=True)
    for organelle in ("chloroplast", "mitochondria"):
        reference = references[organelle]
        analysis_mask = (
            root / "canonical_publication/references/masks/chloroplast_population_sites.bed"
            if organelle == "chloroplast"
            else root / "canonical_publication/references/masks" / args.run_id / "mitochondria_high_confidence_sites.bed"
        )
        analysis_bases = analysis_mask_length(read_mask_intervals(analysis_mask, organelle), len(reference))
        accepted_variants = variant_dir / f"{organelle}.high_confidence_variant_sites.vcf.gz"
        likelihoods = work_dir / f"{organelle}.mpileup_likelihoods.bcf"
        all_sites = work_dir / f"{organelle}.genotype_masked.all_sites.bcf"
        samples_path = metadata_dir / f"{organelle}_samples.tsv"
        variant_state_path = root / "canonical_publication/provenance/runs" / args.run_id / "variants" / f"{organelle}.json"
        if not all(path.exists() for path in (accepted_variants, likelihoods, all_sites, samples_path, variant_state_path)):
            raise RuntimeError(f"Missing canonical variant inputs for {organelle}")
        variant_state = json.loads(variant_state_path.read_text())
        fasta = result_dir / f"{organelle}.callable_alignment.fa"
        summary = result_dir / f"{organelle}.callable_summary.tsv"
        state_path = state_dir / f"{organelle}.json"
        command = (
            f"bcftools query per sample {all_sites.relative_to(root)}; "
            f"derive invariant haploid GQ from homozygous PLs in {likelihoods.relative_to(root)}; "
            f"accept high-confidence variants from {accepted_variants.relative_to(root)}"
        )
        declared = {
            **runtime_provenance(root, {"bcftools": ("bcftools", "--version")}),
            accepted_variants.relative_to(root).as_posix(): sha256_file(accepted_variants),
            likelihoods.relative_to(root).as_posix(): sha256_file(likelihoods),
            all_sites.relative_to(root).as_posix(): sha256_file(all_sites),
            samples_path.relative_to(root).as_posix(): sha256_file(samples_path),
            analysis_mask.relative_to(root).as_posix(): sha256_file(analysis_mask),
            config_path.relative_to(root).as_posix(): sha256_file(config_path),
        }
        fingerprint = build_stage_fingerprint_from_hashes(
            f"consensus:{organelle}",
            declared,
            {"variants": variant_state["fingerprint"]["digest"]},
            [command],
        )
        if args.resume and state_path.exists():
            saved = json.loads(state_path.read_text())
            validate_resume(saved["fingerprint"]["digest"], fingerprint)
            for path, digest in saved["outputs"].items():
                if sha256_file(root / path) != digest:
                    raise RuntimeError(f"Consensus output checksum mismatch: {path}")
            print(f"resume-valid consensus {organelle}")
            continue
        if state_path.exists() or fasta.exists() or summary.exists():
            raise RuntimeError(f"Existing unvalidated consensus output for {organelle}")
        accepted_positions = accepted_variant_positions(accepted_variants)
        records = []
        summary_rows = []
        for sample in read_samples(samples_path):
            sequence = build_callable_sequence(
                reference,
                sample_rows(all_sites, sample),
                accepted_positions,
                invariant_rows=sample_invariant_rows(likelihoods, sample),
                minimum_depth=minimum_depth,
                minimum_genotype_quality=minimum_genotype_quality,
            )
            callable_bases = sum(base in "ACGT" for base in sequence)
            concordance = reference_concordance(reference, sequence)
            records.append((sample, sequence))
            summary_rows.append(
                {
                    "sample_id": sample,
                    "reference_length": len(sequence),
                    "analysis_mask_bases": analysis_bases,
                    "callable_bases": callable_bases,
                    "callable_fraction_of_reference": f"{callable_bases / len(sequence):.8f}",
                    "callable_fraction_of_analysis_mask": f"{callable_bases / analysis_bases:.8f}",
                    "reference_matches_at_callable_bases": concordance.reference_matches,
                    "nonreference_callable_bases": concordance.nonreference_bases,
                    "callable_reference_identity": f"{concordance.identity:.12g}",
                }
            )
            print(f"consensus {organelle} {sample}", flush=True)
        write_fasta(fasta, records)
        with summary.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                delimiter="\t",
                lineterminator="\n",
                fieldnames=[
                    "sample_id",
                    "reference_length",
                    "analysis_mask_bases",
                    "callable_bases",
                    "callable_fraction_of_reference",
                    "callable_fraction_of_analysis_mask",
                    "reference_matches_at_callable_bases",
                    "nonreference_callable_bases",
                    "callable_reference_identity",
                ],
            )
            writer.writeheader()
            writer.writerows(summary_rows)
        outputs = {path.relative_to(root).as_posix(): sha256_file(path) for path in (fasta, summary)}
        state_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "organelle": organelle,
                    "sample_count": len(records),
                    "analysis_mask_bases": analysis_bases,
                    "accepted_high_confidence_variant_site_count": len(accepted_positions),
                    "fingerprint": asdict(fingerprint),
                    "outputs": outputs,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(f"built {organelle} callable alignment for {len(records)} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
