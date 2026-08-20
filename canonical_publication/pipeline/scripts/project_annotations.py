#!/usr/bin/env python3
"""Regenerate draft annotations by BLAST-projecting accessioned source features."""

from __future__ import annotations

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from organelle_pipeline.annotations import gff_phase, projected_interval
from organelle_pipeline.reference_evidence import read_blast_hits
from organelle_pipeline.references import write_fasta

FEATURE_TYPES = {"gene", "CDS", "tRNA", "rRNA"}
BLAST_FORMAT = "6 qseqid sseqid pident length qstart qend sstart send bitscore qlen slen"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, required=True)
    return parser.parse_args()


def clean(value: str) -> str:
    return value.replace(";", ",").replace("=", ":").replace("\t", " ")


def main() -> int:
    args = parse_args()
    root = args.repository_root.resolve()
    annotation_dir = root / "canonical_publication/references/annotations"
    evidence_dir = root / "canonical_publication/references/evidence/annotation_projection"
    annotation_dir.mkdir(parents=True, exist_ok=True)
    evidence_dir.mkdir(parents=True, exist_ok=True)
    sources = {
        "chloroplast": root / "source_data/reference_candidates/external/NC_085682.1.gb",
        "mitochondria": root / "source_data/reference_candidates/external/PV256627.1.gb",
    }
    summaries = []
    for organelle, genbank in sources.items():
        record = SeqIO.read(genbank, "genbank")
        selected = root / f"canonical_publication/references/selected/{organelle}.fa"
        feature_metadata = {}
        query_records = []
        for index, feature in enumerate(record.features, 1):
            if feature.type not in FEATURE_TYPES:
                continue
            sequence = str(feature.extract(record.seq)).upper()
            if len(sequence) < 20 or set(sequence) - set("ACGTN"):
                continue
            query_id = f"feature_{index:05d}"
            qualifiers = feature.qualifiers
            feature_metadata[query_id] = {
                "type": feature.type,
                "gene": (qualifiers.get("gene") or qualifiers.get("locus_tag") or [query_id])[0],
                "product": (qualifiers.get("product") or [""])[0],
                "source_location": str(feature.location),
                "length": len(sequence),
                "codon_start": (qualifiers.get("codon_start") or [None])[0],
            }
            query_records.append((query_id, sequence))
        query_fasta = evidence_dir / f"{organelle}.external_features.fa"
        blast_output = evidence_dir / f"{organelle}.feature_projection.blastn.tsv"
        write_fasta(query_fasta, query_records)
        subprocess.run(
            [
                "blastn",
                "-query",
                str(query_fasta),
                "-subject",
                str(selected),
                "-task",
                "blastn",
                "-dust",
                "no",
                "-evalue",
                "1e-10",
                "-max_target_seqs",
                "5",
                "-max_hsps",
                "1",
                "-outfmt",
                BLAST_FORMAT,
                "-out",
                str(blast_output),
            ],
            check=True,
        )
        grouped = defaultdict(list)
        for hit in read_blast_hits(blast_output.read_text().splitlines()):
            grouped[hit.query].append(hit)
        accepted = []
        rejected = []
        for query_id, metadata in feature_metadata.items():
            hits = grouped.get(query_id, [])
            if not hits:
                rejected.append((query_id, "no_hit"))
                continue
            hit = max(hits, key=lambda value: (value.bitscore, value.alignment_length))
            coverage = hit.alignment_length / hit.query_length
            if coverage < 0.80 or hit.identity_percent < 90.0:
                rejected.append(
                    (
                        query_id,
                        f"below_threshold:coverage={coverage:.4f};identity={hit.identity_percent:.3f}",
                    )
                )
                continue
            start, end, strand = projected_interval(hit)
            accepted.append((query_id, metadata, hit, coverage, start, end, strand))
        gff = annotation_dir / f"{organelle}.projected.gff3"
        table = annotation_dir / f"{organelle}.projected.tsv"
        with gff.open("w") as handle:
            handle.write("##gff-version 3\n")
            handle.write(
                f"##provenance external={genbank.relative_to(root)}; method=feature_blast_projection; "
                "minimum_query_coverage=0.80; minimum_identity=90.0; status=draft\n"
            )
            for query_id, metadata, hit, coverage, start, end, strand in sorted(accepted, key=lambda value: (value[4], value[5], value[0])):
                attributes = (
                    f"ID={query_id};gene={clean(metadata['gene'])};"
                    f"product={clean(metadata['product'])};external_accession={record.id};"
                    f"query_coverage={coverage:.6f};identity={hit.identity_percent:.6f};"
                    "status=draft_projection"
                )
                handle.write(
                    f"{organelle}\tcanonical_projection\t{metadata['type']}\t{start + 1}\t{end}\t"
                    f"{hit.bitscore:.3f}\t{strand}\t{gff_phase(metadata['type'], metadata['codon_start'])}\t"
                    f"{attributes}\n"
                )
        with table.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "feature_id",
                    "feature_type",
                    "gene",
                    "product",
                    "start_1based",
                    "end_1based",
                    "strand",
                    "query_coverage",
                    "identity_percent",
                    "external_accession",
                    "status",
                ]
            )
            for query_id, metadata, hit, coverage, start, end, strand in sorted(accepted, key=lambda value: (value[4], value[5], value[0])):
                writer.writerow(
                    [
                        query_id,
                        metadata["type"],
                        metadata["gene"],
                        metadata["product"],
                        start + 1,
                        end,
                        strand,
                        f"{coverage:.6f}",
                        f"{hit.identity_percent:.6f}",
                        record.id,
                        "draft_projection",
                    ]
                )
        rejection_path = evidence_dir / f"{organelle}.unmapped_features.tsv"
        with rejection_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["feature_id", "reason"])
            writer.writerows(rejected)
        summaries.append(
            {
                "organelle": organelle,
                "external_accession": record.id,
                "candidate_features": len(feature_metadata),
                "mapped_features": len(accepted),
                "unmapped_features": len(rejected),
                "mapped_fraction": len(accepted) / len(feature_metadata),
                "status": "draft_projection_not_de_novo_annotation",
            }
        )
    summary_path = evidence_dir / "annotation_projection_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            lineterminator="\n",
            fieldnames=list(summaries[0]),
        )
        writer.writeheader()
        writer.writerows(summaries)
    print("projected draft annotations from independently sourced GenBank records")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
