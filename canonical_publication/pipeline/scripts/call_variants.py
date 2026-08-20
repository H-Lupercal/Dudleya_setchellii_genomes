#!/usr/bin/env python3
"""Call all-site haploid genotypes and canonical SNP sets per organelle."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import tomllib
from dataclasses import asdict
from pathlib import Path

from organelle_pipeline.paths import repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
)
from organelle_pipeline.variants import (
    build_all_sites_call_command,
    build_primary_filter_commands,
)


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


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    reference = root / "canonical_publication/references/selected/organelle_combined.fa"
    qc_state_path = root / "canonical_publication/provenance/runs" / args.run_id / "qc.json"
    if not qc_state_path.exists():
        raise RuntimeError("Variant calling requires completed canonical QC")
    qc_state = json.loads(qc_state_path.read_text())
    work = root / "canonical_publication/work" / args.run_id / "variants"
    results = root / "canonical_publication/results/variants" / args.run_id
    metadata = root / "canonical_publication/metadata/qc" / args.run_id
    state_dir = root / "canonical_publication/provenance/runs" / args.run_id / "variants"
    log_dir = root / "canonical_publication/provenance/runs" / args.run_id / "logs/variants"
    for directory in (work, results, state_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)
    organelles = {
        "chloroplast": root / "canonical_publication/references/masks/chloroplast_population_sites.bed",
        "mitochondria": root / "canonical_publication/references/masks" / args.run_id / "mitochondria_high_confidence_sites.bed",
    }

    def relative(path: Path) -> Path:
        return path.relative_to(root)

    for organelle, mask in organelles.items():
        samples = read_samples(metadata / f"{organelle}_samples.tsv")
        bam_list = work / f"{organelle}.bam.list"
        bam_list.write_text("".join(f"canonical_publication/work/{args.run_id}/mapping/{sample}.organelle.bam\n" for sample in samples))
        likelihoods = work / f"{organelle}.mpileup_likelihoods.bcf"
        all_sites = work / f"{organelle}.all_sites.bcf"
        masked = work / f"{organelle}.genotype_masked.all_sites.bcf"
        high_confidence = results / f"{organelle}.high_confidence_variant_sites.vcf.gz"
        primary = results / f"{organelle}.primary.vcf.gz"
        mac2 = results / f"{organelle}.mac2_ordination.vcf.gz"
        commands = (
            build_all_sites_call_command(
                relative(reference),
                relative(bam_list),
                relative(mask),
                relative(all_sites),
                likelihood_bcf=relative(likelihoods),
                minimum_mapping_quality=int(config["mapping"]["minimum_mapping_quality"]),
                minimum_base_quality=int(config["mapping"]["minimum_base_quality"]),
                maximum_per_file_depth=int(config["variants"]["maximum_per_file_pileup_depth"]),
                ploidy=int(config["variants"]["ploidy"]),
            ),
            *build_primary_filter_commands(
                relative(reference),
                relative(all_sites),
                relative(masked),
                relative(high_confidence),
                relative(primary),
                relative(mac2),
                minimum_depth=int(config["variants"]["minimum_depth"]),
                minimum_genotype_quality=int(config["variants"]["minimum_genotype_quality"]),
                minimum_site_quality=int(config["variants"]["minimum_site_quality"]),
                maximum_missing_fraction=float(config["variants"]["maximum_missing_fraction"]),
                primary_minimum_mac=int(config["variants"]["primary_minimum_minor_allele_count"]),
                ordination_minimum_mac=int(config["variants"]["ordination_minimum_minor_allele_count"]),
            ),
        )
        declared = {
            **runtime_provenance(root, {"bcftools": ("bcftools", "--version")}),
            relative(reference).as_posix(): sha256_file(reference),
            relative(mask).as_posix(): sha256_file(mask),
            relative(bam_list).as_posix(): sha256_file(bam_list),
            relative(config_path).as_posix(): sha256_file(config_path),
        }
        fingerprint = build_stage_fingerprint_from_hashes(
            f"variants:{organelle}",
            declared,
            {"qc": qc_state["fingerprint"]["digest"]},
            commands,
        )
        state_path = state_dir / f"{organelle}.json"
        if args.resume and state_path.exists():
            saved = json.loads(state_path.read_text())
            validate_resume(saved["fingerprint"]["digest"], fingerprint)
            for path, digest in saved["outputs"].items():
                if sha256_file(root / path) != digest:
                    raise RuntimeError(f"Variant output checksum mismatch: {path}")
            print(f"resume-valid variants {organelle}")
            continue
        if state_path.exists() or any(path.exists() for path in (likelihoods, all_sites, masked, high_confidence, primary, mac2)):
            raise RuntimeError(f"Existing unvalidated variant output for {organelle}")
        log_path = log_dir / f"{organelle}.log"
        with log_path.open("w") as log:
            for command in commands:
                log.write(f"COMMAND\t{command}\n")
                log.flush()
                subprocess.run(
                    ["bash", "-o", "pipefail", "-c", command],
                    cwd=root,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    env=os.environ.copy(),
                    check=True,
                )
        outputs = {
            relative(path).as_posix(): sha256_file(path)
            for path in (
                likelihoods,
                all_sites,
                masked,
                high_confidence,
                Path(f"{high_confidence}.csi"),
                primary,
                Path(f"{primary}.csi"),
                mac2,
                Path(f"{mac2}.csi"),
            )
        }
        state_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "organelle": organelle,
                    "sample_count": len(samples),
                    "fingerprint": asdict(fingerprint),
                    "outputs": outputs,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(f"called {organelle} variants for {len(samples)} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
