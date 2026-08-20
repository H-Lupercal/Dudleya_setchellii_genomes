#!/usr/bin/env python3
"""Run explicitly supplementary pseudo-diploid ADMIXTURE replicates."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import tomllib
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

from organelle_pipeline.admixture import pseudo_diploid_alleles, validate_q_matrix
from organelle_pipeline.analysis import select_best_k
from organelle_pipeline.logs import portable_command_log
from organelle_pipeline.paths import repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)

CV_PATTERN = re.compile(r"CV error \(K=(\d+)\):\s*([0-9.eE+-]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--threads-per-job", type=int, default=2)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def vcf_matrix(vcf: Path) -> tuple[list[str], list[tuple[int, str, str, list[str]]]]:
    samples = subprocess.run(
        ["bcftools", "query", "-l", str(vcf)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.splitlines()
    completed = subprocess.run(
        ["bcftools", "query", "-f", "%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)],
        capture_output=True,
        text=True,
        check=True,
    )
    variants = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        variants.append((int(fields[0]), fields[1], fields[2], fields[3:]))
    return samples, variants


def write_ped_map(prefix: Path, samples: list[str], variants: list[tuple[int, str, str, list[str]]]) -> None:
    with Path(f"{prefix}.map").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=" ", lineterminator="\n")
        for index, (position, _, _, _) in enumerate(variants, 1):
            writer.writerow([1, f"organelle_{position}_{index}", 0, position])
    with Path(f"{prefix}.ped").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=" ", lineterminator="\n")
        for sample_index, sample in enumerate(samples):
            row = [sample, sample, 0, 0, 0, -9]
            for _, ref, alt, genotypes in variants:
                row.extend(pseudo_diploid_alleles(genotypes[sample_index], ref, alt))
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    minimum_k = int(config["admixture"]["minimum_k"])
    maximum_k = int(config["admixture"]["maximum_k"])
    replicate_count = int(config["admixture"]["replicates"])
    seed_base = int(config["admixture"]["seed_base"])
    variant_dir = root / "canonical_publication/results/variants" / args.run_id
    work_dir = root / "canonical_publication/work" / args.run_id / "admixture"
    result_dir = root / "canonical_publication/results/supplement" / args.run_id / "admixture"
    state_dir = root / "canonical_publication/provenance/runs" / args.run_id / "admixture"
    log_dir = root / "canonical_publication/provenance/runs" / args.run_id / "logs/admixture"
    for directory in (work_dir, result_dir, state_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)
    for organelle in ("chloroplast", "mitochondria"):
        vcf = variant_dir / f"{organelle}.mac2_ordination.vcf.gz"
        variant_state = json.loads(
            (root / "canonical_publication/provenance/runs" / args.run_id / "variants" / f"{organelle}.json").read_text()
        )
        commands = [
            f"admixture --cv=10 --seed={seed_base + k * 100 + replicate} -j{args.threads_per_job} input.bed {k}"
            for k in range(minimum_k, maximum_k + 1)
            for replicate in range(1, replicate_count + 1)
        ]
        fingerprint = build_stage_fingerprint_from_hashes(
            f"admixture_supplementary:{organelle}",
            {
                **runtime_provenance(
                    root,
                    {
                        "admixture": ("admixture", "--version"),
                        "bcftools": ("bcftools", "--version"),
                        "plink": ("plink", "--version"),
                    },
                ),
                config_path.relative_to(root).as_posix(): sha256_file(config_path),
                vcf.relative_to(root).as_posix(): sha256_file(vcf),
            },
            {"variants": variant_state["fingerprint"]["digest"]},
            commands,
        )
        organelle_result = result_dir / organelle
        organelle_work = work_dir / organelle
        organelle_result.mkdir(parents=True, exist_ok=True)
        organelle_work.mkdir(parents=True, exist_ok=True)
        state_path = state_dir / f"{organelle}.json"
        if args.resume and state_path.exists():
            saved = json.loads(state_path.read_text())
            validate_resume(saved["fingerprint"]["digest"], fingerprint)
            validate_saved_outputs(root, saved)
            print(f"resume-valid supplementary ADMIXTURE {organelle}")
            continue
        if state_path.exists() or any(organelle_result.iterdir()) or any(organelle_work.iterdir()):
            raise RuntimeError(f"Existing unvalidated supplementary ADMIXTURE output for {organelle}")
        samples, variants = vcf_matrix(vcf)
        if not variants:
            raise RuntimeError(f"No MAC>=2 markers available for ADMIXTURE: {vcf}")
        sample_order_path = organelle_result / "sample_order.tsv"
        with sample_order_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["q_row_1based", "sample_id"])
            writer.writerows((index, sample) for index, sample in enumerate(samples, 1))
        input_prefix = organelle_work / "input"
        write_ped_map(input_prefix, samples, variants)
        plink_log_path = log_dir / f"{organelle}.plink.log"
        plink_completed = subprocess.run(
            [
                "plink",
                "--file",
                str(input_prefix),
                "--make-bed",
                "--allow-no-sex",
                "--out",
                str(input_prefix),
            ],
            cwd=root,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=True,
        )
        plink_log_path.write_text(portable_command_log(plink_completed.stdout, root))
        fam_sample_order = [line.split()[1] for line in Path(f"{input_prefix}.fam").read_text().splitlines() if line.strip()]
        if fam_sample_order != samples:
            raise RuntimeError(f"PLINK changed sample order for {organelle}; refusing to label Q rows")

        def run_replicate(
            k: int,
            replicate: int,
            sample_count: int,
            organelle: str = organelle,
            organelle_work: Path = organelle_work,
            input_prefix: Path = input_prefix,
            organelle_result: Path = organelle_result,
        ) -> dict[str, object]:
            seed = seed_base + k * 100 + replicate
            rep_dir = organelle_work / f"K{k}" / f"replicate_{replicate}"
            rep_dir.mkdir(parents=True, exist_ok=True)
            for extension in ("bed", "bim", "fam"):
                source = Path(f"{input_prefix}.{extension}")
                destination = rep_dir / f"input.{extension}"
                if not destination.exists():
                    os.link(source, destination)
            command = [
                "admixture",
                "--cv=10",
                f"--seed={seed}",
                f"-j{args.threads_per_job}",
                "input.bed",
                str(k),
            ]
            completed = subprocess.run(
                command,
                cwd=rep_dir,
                capture_output=True,
                text=True,
                check=True,
            )
            text = completed.stdout + "\n" + completed.stderr
            (log_dir / f"{organelle}.K{k}.replicate_{replicate}.log").write_text(portable_command_log(text, root))
            match = CV_PATTERN.search(text)
            if match is None:
                raise RuntimeError(f"No CV error found for {organelle} K{k} replicate {replicate}")
            q_source = rep_dir / f"input.{k}.Q"
            validate_q_matrix(q_source.read_text().splitlines(), sample_count=sample_count, k=k)
            q_destination = organelle_result / f"K{k}.replicate_{replicate}.Q.tsv"
            shutil.copyfile(q_source, q_destination)
            return {
                "organelle": organelle,
                "k": k,
                "replicate": replicate,
                "seed": seed,
                "cv_error": float(match.group(2)),
                "q_path": q_destination.relative_to(root).as_posix(),
            }

        rows = []
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = [
                executor.submit(run_replicate, k, replicate, len(samples))
                for k in range(minimum_k, maximum_k + 1)
                for replicate in range(1, replicate_count + 1)
            ]
            for future in as_completed(futures):
                rows.append(future.result())
                print(
                    f"admixture {organelle} K{rows[-1]['k']} replicate {rows[-1]['replicate']}",
                    flush=True,
                )
        rows.sort(key=lambda row: (int(row["k"]), int(row["replicate"])))
        replicate_path = organelle_result / "replicate_cv.tsv"
        with replicate_path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                delimiter="\t",
                lineterminator="\n",
                fieldnames=["organelle", "k", "replicate", "seed", "cv_error", "q_path"],
            )
            writer.writeheader()
            writer.writerows(rows)
        means = {k: sum(float(row["cv_error"]) for row in rows if row["k"] == k) / replicate_count for k in range(minimum_k, maximum_k + 1)}
        choice = select_best_k(means, minimum_k, maximum_k)
        selected_replicate = min(
            (row for row in rows if row["k"] == choice.k),
            key=lambda row: (float(row["cv_error"]), int(row["replicate"])),
        )
        selected_q_path = organelle_result / f"selected_K{choice.k}.best_cv.Q.tsv"
        shutil.copyfile(root / str(selected_replicate["q_path"]), selected_q_path)
        selected_solution_path = organelle_result / "selected_solution.tsv"
        with selected_solution_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["selected_k", "replicate", "seed", "cv_error", "boundary_optimum", "q_path"])
            writer.writerow(
                [
                    choice.k,
                    selected_replicate["replicate"],
                    selected_replicate["seed"],
                    selected_replicate["cv_error"],
                    "yes" if choice.is_boundary else "no",
                    selected_q_path.relative_to(root).as_posix(),
                ]
            )
        summary_path = organelle_result / "k_summary.tsv"
        with summary_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["k", "mean_cv_error", "selected", "boundary_optimum"])
            for k in range(minimum_k, maximum_k + 1):
                writer.writerow(
                    [
                        k,
                        f"{means[k]:.12g}",
                        "yes" if k == choice.k else "no",
                        "yes" if k == choice.k and choice.is_boundary else "no",
                    ]
                )
        limitation_path = organelle_result / "INTERPRETATION_LIMITATIONS.md"
        limitation_path.write_text(
            "# Supplementary ADMIXTURE limitations\n\n"
            "These maternally inherited haploid organelle SNPs were duplicated into "
            "homozygous pseudo-diploid genotypes solely for software compatibility. "
            "Markers are physically linked, were not LD-pruned, and violate the "
            "autosomal/unlinked assumptions of ADMIXTURE. Results are descriptive "
            "haplotype-clustering sensitivity analyses, not ancestry proportions.\n"
        )
        q_outputs = [root / str(row["q_path"]) for row in rows]
        log_outputs = [
            plink_log_path,
            *[
                log_dir / f"{organelle}.K{k}.replicate_{replicate}.log"
                for k in range(minimum_k, maximum_k + 1)
                for replicate in range(1, replicate_count + 1)
            ],
        ]
        outputs = {
            path.relative_to(root).as_posix(): sha256_file(path)
            for path in (
                replicate_path,
                summary_path,
                limitation_path,
                sample_order_path,
                selected_q_path,
                selected_solution_path,
                *q_outputs,
                *log_outputs,
            )
        }
        state_path.write_text(
            json.dumps(
                {
                    "status": "complete",
                    "role": "supplementary",
                    "organelle": organelle,
                    "selected_k": choice.k,
                    "boundary_optimum": choice.is_boundary,
                    "validated_q_matrix_count": len(q_outputs),
                    "fingerprint": asdict(fingerprint),
                    "outputs": outputs,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(f"completed supplementary ADMIXTURE {organelle}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
