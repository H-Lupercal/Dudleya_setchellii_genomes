#!/usr/bin/env python3
"""Reconstruct and independently verify canonical organelle references."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tomllib
from dataclasses import asdict
from pathlib import Path

from Bio import SeqIO
from organelle_pipeline.inventory import ACCEPTABLE_SOURCE_VALIDATION_STATUSES
from organelle_pipeline.paths import assert_canonical_path, repository_relative, validate_run_id
from organelle_pipeline.provenance import (
    build_stage_fingerprint_from_hashes,
    runtime_provenance,
    sha256_file,
    validate_resume,
    validate_saved_outputs,
)
from organelle_pipeline.reference_evidence import (
    read_blast_hits,
    select_inverted_repeat_pair,
    select_terminal_direct_repeat,
    self_repeat_intervals,
    summarize_blast_hits,
    summarize_query_orientation,
)
from organelle_pipeline.references import (
    normalize_chloroplast,
    read_single_fasta,
    write_combined_reference,
    write_fasta,
)

BLAST_FORMAT = "6 qseqid sseqid pident length qstart qend sstart send bitscore qlen slen"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def run_blast(query: Path, subject: Path, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "blastn",
            "-query",
            str(query),
            "-subject",
            str(subject),
            "-task",
            "blastn",
            "-dust",
            "no",
            "-evalue",
            "1e-20",
            "-outfmt",
            BLAST_FORMAT,
            "-out",
            str(output),
        ],
        check=True,
    )


def write_bed(path: Path, rows: list[tuple[str, int, int, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


def main() -> int:
    args = parse_args()
    validate_run_id(args.run_id)
    root = args.repository_root.resolve()
    config_path = root / repository_relative(args.config, root)
    config = tomllib.loads(config_path.read_text())
    reference_config = config["references"]

    def source(name: str) -> Path:
        return assert_canonical_path(root / reference_config[name], root)

    cp_candidate = source("chloroplast_candidate")
    mt_candidate = source("mitochondria_candidate")
    cp_external = source("chloroplast_external")
    mt_external = source("mitochondria_external")
    cp_genbank = root / "source_data/reference_candidates/external/NC_085682.1.gb"
    mt_genbank = root / "source_data/reference_candidates/external/PV256627.1.gb"
    selected = root / "canonical_publication/references/selected"
    evidence = root / "canonical_publication/references/evidence"
    masks = root / "canonical_publication/references/masks"
    annotations = root / "canonical_publication/references/annotations"
    state_path = root / "canonical_publication/provenance/runs" / args.run_id / "references.json"
    source_state_path = root / "canonical_publication/provenance/runs" / args.run_id / "source_validation.json"
    if not source_state_path.is_file():
        raise RuntimeError("Reference preparation requires completed immutable-source validation")
    source_state = json.loads(source_state_path.read_text())
    if source_state.get("status") not in ACCEPTABLE_SOURCE_VALIDATION_STATUSES:
        raise RuntimeError("Immutable-source validation did not pass")
    validate_saved_outputs(root, source_state)
    selected.mkdir(parents=True, exist_ok=True)
    evidence.mkdir(parents=True, exist_ok=True)
    masks.mkdir(parents=True, exist_ok=True)
    state_path.parent.mkdir(parents=True, exist_ok=True)

    declared_inputs = {
        **runtime_provenance(
            root,
            {
                "biopython": ("python", "-c", "import Bio; print(Bio.__version__)"),
                "blastn": ("blastn", "-version"),
                "bwa": ("bwa",),
                "samtools": ("samtools", "--version"),
            },
        ),
        config_path.relative_to(root).as_posix(): sha256_file(config_path),
        source_state_path.relative_to(root).as_posix(): sha256_file(source_state_path),
        cp_candidate.relative_to(root).as_posix(): sha256_file(cp_candidate),
        mt_candidate.relative_to(root).as_posix(): sha256_file(mt_candidate),
        cp_external.relative_to(root).as_posix(): sha256_file(cp_external),
        mt_external.relative_to(root).as_posix(): sha256_file(mt_external),
    }
    for genbank in (cp_genbank, mt_genbank):
        declared_inputs[genbank.relative_to(root).as_posix()] = sha256_file(genbank)
    fingerprint = build_stage_fingerprint_from_hashes(
        "references",
        declared_inputs,
        {"immutable_sources": source_state["fingerprint"]["digest"]},
        [
            "blastn candidate/external and self comparisons; validate chloroplast terminal repeat and IRs; mask mitochondrial repeats",
            "project external GenBank annotations by BLAST; bwa index; samtools faidx",
        ],
    )
    if args.resume and state_path.exists():
        saved = json.loads(state_path.read_text())
        validate_resume(saved["fingerprint"]["digest"], fingerprint)
        for path, digest in saved["outputs"].items():
            if sha256_file(root / path) != digest:
                raise RuntimeError(f"Reference output checksum mismatch: {path}")
        print("resume-valid references")
        return 0
    if state_path.exists():
        raise RuntimeError("Reference state already exists; use --resume or a new run ID")

    with (evidence / "external_accession_sequence_consistency.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "organelle",
                "fasta_accession",
                "genbank_accession",
                "fasta_length",
                "genbank_length",
                "sequence_identical",
            ]
        )
        for organelle, fasta_path, genbank_path in (
            ("chloroplast", cp_external, cp_genbank),
            ("mitochondria", mt_external, mt_genbank),
        ):
            fasta_accession, fasta_sequence = read_single_fasta(fasta_path)[0]
            genbank_record = SeqIO.read(genbank_path, "genbank")
            genbank_sequence = str(genbank_record.seq).upper()
            identical = fasta_accession == genbank_record.id and fasta_sequence == genbank_sequence
            writer.writerow(
                [
                    organelle,
                    fasta_accession,
                    genbank_record.id,
                    len(fasta_sequence),
                    len(genbank_sequence),
                    "yes" if identical else "no",
                ]
            )
            if not identical:
                raise RuntimeError(
                    f"External FASTA/GenBank accession or sequence mismatch for {organelle}: {fasta_accession} vs {genbank_record.id}"
                )

    cp_raw = read_single_fasta(cp_candidate)[0][1]
    mt_sequence = read_single_fasta(mt_candidate)[0][1]
    deduplicated_length = int(reference_config["chloroplast_deduplicated_length"])
    raw_self_path = evidence / "chloroplast_raw_self.blastn.tsv"
    run_blast(cp_candidate, cp_candidate, raw_self_path)
    raw_self_hits = read_blast_hits(raw_self_path.read_text().splitlines())
    terminal_repeat = select_terminal_direct_repeat(
        raw_self_hits,
        retained_length=deduplicated_length,
        sequence_length=len(cp_raw),
        minimum_identity=float(reference_config["chloroplast_terminal_duplicate_minimum_identity_percent"]),
    )
    cp_sequence, removed = normalize_chloroplast(cp_raw, deduplicated_length)
    cp_path = selected / "chloroplast.fa"
    mt_path = selected / "mitochondria.fa"
    combined_path = selected / "organelle_combined.fa"
    write_fasta(cp_path, [("chloroplast", cp_sequence)])
    write_fasta(mt_path, [("mitochondria", mt_sequence)])
    write_combined_reference(combined_path, cp_sequence, mt_sequence)
    write_fasta(
        evidence / "chloroplast_removed_terminal_sequence.fa",
        [(f"removed_after_{deduplicated_length}", removed)],
    )
    with (evidence / "chloroplast_terminal_repeat_validation.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["copy", "start_0based", "end_0based_exclusive", "interpretation"])
        writer.writerow(["start_copy", *terminal_repeat[0], "self-BLAST-supported terminal direct repeat"])
        writer.writerow(["end_copy", *terminal_repeat[1], "trimmed redundant terminal copy"])

    comparisons = {
        "chloroplast_raw_vs_NC_085682.1": (cp_candidate, cp_external),
        "chloroplast_selected_vs_NC_085682.1": (cp_path, cp_external),
        "mitochondria_selected_vs_PV256627.1": (mt_path, mt_external),
        "chloroplast_selected_self": (cp_path, cp_path),
        "mitochondria_selected_self": (mt_path, mt_path),
    }
    summaries = {"chloroplast_raw_self": summarize_blast_hits(raw_self_hits)}
    blast_paths = {"chloroplast_raw_self": raw_self_path}
    for label, (query, subject) in comparisons.items():
        blast_path = evidence / f"{label}.blastn.tsv"
        run_blast(query, subject, blast_path)
        hits = read_blast_hits(blast_path.read_text().splitlines())
        summaries[label] = summarize_blast_hits(hits)
        blast_paths[label] = blast_path

    with (evidence / "reference_orientation_and_boundary_checks.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "organelle",
                "dominant_external_hsp_orientation",
                "dominant_hsp_identity_percent",
                "dominant_hsp_aligned_bp",
                "same_orientation_assigned_query_bp",
                "reverse_complement_assigned_query_bp",
                "same_orientation_fraction_of_covered_query",
                "reverse_complement_fraction_of_covered_query",
                "assembly_boundary_evidence",
            ]
        )
        for organelle, label in (
            ("chloroplast", "chloroplast_selected_vs_NC_085682.1"),
            ("mitochondria", "mitochondria_selected_vs_PV256627.1"),
        ):
            hits = read_blast_hits(blast_paths[label].read_text().splitlines())
            dominant = max(hits, key=lambda hit: (hit.bitscore, hit.alignment_length))
            same_orientation = (dominant.query_end - dominant.query_start) * (dominant.subject_end - dominant.subject_start) > 0
            orientation = summarize_query_orientation(hits)
            boundary = (
                "terminal direct-repeat trim validated by raw self-BLAST; read-backed check pending QC"
                if organelle == "chloroplast"
                else "candidate retained intact; read-backed high-confidence boundary check pending QC"
            )
            writer.writerow(
                [
                    organelle,
                    "same" if same_orientation else "reverse_complement",
                    dominant.identity_percent,
                    dominant.alignment_length,
                    orientation.same_orientation_query_bp,
                    orientation.reverse_complement_query_bp,
                    f"{orientation.same_orientation_fraction_of_covered:.8f}",
                    f"{orientation.reverse_complement_fraction_of_covered:.8f}",
                    boundary,
                ]
            )

    cp_summary = summaries["chloroplast_selected_vs_NC_085682.1"]
    mt_summary = summaries["mitochondria_selected_vs_PV256627.1"]
    if cp_summary.query_coverage < float(reference_config["minimum_chloroplast_query_coverage"]):
        raise RuntimeError("Selected chloroplast failed query-coverage validation")
    if mt_summary.query_coverage < float(reference_config["minimum_mitochondria_query_coverage"]):
        raise RuntimeError("Selected mitochondrion failed query-coverage validation")
    minimum_identity = float(reference_config["minimum_identity_percent"])
    if cp_summary.weighted_identity_percent < minimum_identity:
        raise RuntimeError("Selected chloroplast failed identity validation")
    if mt_summary.weighted_identity_percent < minimum_identity:
        raise RuntimeError("Selected mitochondrion failed identity validation")

    ir_first, ir_second = select_inverted_repeat_pair(read_blast_hits(blast_paths["chloroplast_selected_self"].read_text().splitlines()))
    write_bed(
        masks / "chloroplast_ir_copies.bed",
        [
            ("chloroplast", ir_first[0], ir_first[1], "IR_copy_1"),
            ("chloroplast", ir_second[0], ir_second[1], "IR_copy_2"),
        ],
    )
    write_bed(
        masks / "chloroplast_duplicate_ir_mask.bed",
        [("chloroplast", ir_second[0], ir_second[1], "duplicate_IR_copy")],
    )
    unique_eligibility_intervals = []
    cursor = 0
    for start, end in (ir_first, ir_second):
        if cursor < start:
            unique_eligibility_intervals.append(("chloroplast", cursor, start, "unique_eligibility_sites"))
        cursor = max(cursor, end)
    if cursor < len(cp_sequence):
        unique_eligibility_intervals.append(("chloroplast", cursor, len(cp_sequence), "unique_eligibility_sites"))
    write_bed(
        masks / "chloroplast_unique_eligibility_sites.bed",
        unique_eligibility_intervals,
    )
    population_intervals = []
    if ir_second[0] > 0:
        population_intervals.append(("chloroplast", 0, ir_second[0], "population_sites_before_duplicate_IR"))
    if ir_second[1] < len(cp_sequence):
        population_intervals.append(
            (
                "chloroplast",
                ir_second[1],
                len(cp_sequence),
                "population_sites_after_duplicate_IR",
            )
        )
    write_bed(masks / "chloroplast_population_sites.bed", population_intervals)
    write_bed(
        masks / "chloroplast_full_reference.bed",
        [("chloroplast", 0, len(cp_sequence), "full_reference")],
    )
    write_bed(
        masks / "mitochondria_full_reference.bed",
        [("mitochondria", 0, len(mt_sequence), "full_reference")],
    )
    mt_repeats = self_repeat_intervals(
        read_blast_hits(blast_paths["mitochondria_selected_self"].read_text().splitlines()),
        minimum_length=int(reference_config["mitochondria_repeat_minimum_length"]),
        minimum_identity=float(reference_config["mitochondria_repeat_minimum_identity_percent"]),
    )
    write_bed(
        masks / "mitochondria_repeat_mask.bed",
        [("mitochondria", start, end, f"self_repeat_{index}") for index, (start, end) in enumerate(mt_repeats, 1)],
    )
    mt_repeat_bases = sum(end - start for start, end in mt_repeats)
    with (evidence / "mitochondria_repeat_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "reference_length",
                "repeat_intervals",
                "repeat_masked_bases",
                "repeat_masked_fraction",
                "remaining_unique_bases_before_read_qc",
                "interpretation",
            ]
        )
        writer.writerow(
            [
                len(mt_sequence),
                len(mt_repeats),
                mt_repeat_bases,
                f"{mt_repeat_bases / len(mt_sequence):.8f}",
                len(mt_sequence) - mt_repeat_bases,
                "both copies of self-alignments meeting configured length/identity thresholds are excluded; candidate retained intact",
            ]
        )

    with (evidence / "reference_similarity_summary.tsv").open("w", newline="") as handle:
        fields = ["comparison", *asdict(cp_summary).keys()]
        writer = csv.DictWriter(handle, delimiter="\t", lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        for label, summary in summaries.items():
            writer.writerow({"comparison": label, **asdict(summary)})
    with (evidence / "reference_selection.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["organelle", "raw_length", "selected_length", "removed_length"])
        writer.writerow(["chloroplast", len(cp_raw), len(cp_sequence), len(removed)])
        writer.writerow(["mitochondria", len(mt_sequence), len(mt_sequence), 0])

    subprocess.run(["bwa", "index", str(combined_path)], check=True)
    subprocess.run(["samtools", "faidx", str(combined_path)], check=True)
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).with_name("project_annotations.py")),
            "--repository-root",
            str(root),
        ],
        check=True,
    )
    output_files = sorted(
        {
            *(path for base in (selected, annotations) for path in base.rglob("*") if path.is_file()),
            *(path for path in masks.iterdir() if path.is_file()),
            *(path for path in evidence.iterdir() if path.is_file()),
            *(path for path in (evidence / "annotation_projection").rglob("*") if path.is_file()),
        }
    )
    state_path.write_text(
        json.dumps(
            {
                "status": "complete",
                "fingerprint": asdict(fingerprint),
                "outputs": {path.relative_to(root).as_posix(): sha256_file(path) for path in output_files},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"selected chloroplast {len(cp_sequence)} bp; mitochondria {len(mt_sequence)} bp")
    print(f"chloroplast IR intervals: {ir_first}, {ir_second}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
