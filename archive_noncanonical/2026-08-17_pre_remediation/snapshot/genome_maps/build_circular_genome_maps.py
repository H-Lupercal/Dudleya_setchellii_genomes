#!/usr/bin/env python3
"""Build circular cpDNA and mtDNA genome maps from the reference GFF3 files."""

from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("genome_maps") / ".matplotlib"))

import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from matplotlib.patches import Patch
from pycirclize import Circos
from pycirclize.parser import Gff


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "genome_maps"

INPUTS = {
    "cpDNA": {
        "title": "Dudleya setchellii chloroplast genome",
        "fasta": ROOT
        / "dudleya_organelle_reference_verification"
        / "references"
        / "chloroplast.normalized.fa",
        "gff": ROOT
        / "dudleya_organelle_reference_verification"
        / "annotations"
        / "chloroplast.gff3",
        "snp_sites": ROOT
        / "full_pipeline_run"
        / "results"
        / "10_snp_alignment"
        / "cpDNA.primary.snp_sites.tsv",
        "seqid": "chloroplast",
        "tick_interval": 25_000,
    },
    "mtDNA": {
        "title": "Dudleya setchellii mitochondrial genome",
        "fasta": ROOT
        / "dudleya_organelle_reference_verification"
        / "references"
        / "mitochondria.fa",
        "gff": ROOT
        / "dudleya_organelle_reference_verification"
        / "annotations"
        / "mitochondria.gff3",
        "snp_sites": ROOT
        / "full_pipeline_run"
        / "results"
        / "10_snp_alignment"
        / "mtDNA.primary.snp_sites.tsv",
        "seqid": "mitochondria",
        "tick_interval": 50_000,
    },
}


def read_single_fasta(path: Path):
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError(f"Expected one FASTA record in {path}, found {len(records)}")
    return records[0]


def feature_name(feature: SeqFeature) -> str:
    return feature.qualifiers.get("Name", feature.qualifiers.get("ID", ["unknown"]))[0]


def confidence(feature: SeqFeature) -> str:
    return feature.qualifiers.get("confidence", ["unknown"])[0]


def feature_len(feature: SeqFeature) -> int:
    return int(feature.location.end) - int(feature.location.start)


def color_for_cds(feature: SeqFeature) -> str:
    name = feature_name(feature).lower()
    conf = confidence(feature)
    if conf in {"low", "weak"}:
        return "#b8b8b8"
    if name.startswith(("psa", "psb", "pet", "rbc", "lhb")):
        return "#2f80ed"
    if name.startswith(("atp", "ndh")):
        return "#27ae60"
    if name.startswith(("rpl", "rps", "rrn")):
        return "#9b59b6"
    if name.startswith(("rpo", "mat")):
        return "#f2994a"
    if name.startswith(("cox", "cob", "ccm", "nad")):
        return "#d35400"
    return "#34495e"


def color_for_gene(feature: SeqFeature) -> str:
    conf = confidence(feature)
    if conf in {"low", "weak"}:
        return "#b8b8b8"
    return "#4b5563"


def sliding_gc(seq: str, window: int) -> tuple[list[int], list[float], float]:
    seq = seq.upper()
    genome_gc = (seq.count("G") + seq.count("C")) / len(seq)
    positions: list[int] = []
    values: list[float] = []
    step = max(500, window // 5)
    for start in range(0, len(seq), step):
        window_seq = (seq + seq[:window])[start : start + window]
        gc = (window_seq.count("G") + window_seq.count("C")) / len(window_seq)
        positions.append(min(start + window // 2, len(seq) - 1))
        values.append(gc - genome_gc)
    return positions, values, genome_gc


def read_snp_positions(path: Path) -> list[int]:
    positions: list[int] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            positions.append(int(row["position"]))
    return positions


def snp_density(positions: list[int], genome_len: int, bin_size: int) -> tuple[list[int], list[int]]:
    bin_count = genome_len // bin_size + 1
    counts = [0 for _ in range(bin_count)]
    for position in positions:
        if position < 1:
            continue
        index = min((position - 1) // bin_size, bin_count - 1)
        counts[index] += 1
    centers = [min(i * bin_size + bin_size // 2, genome_len - 1) for i in range(bin_count)]
    return centers, counts


def build_map(label: str, config: dict[str, object]) -> dict[str, object]:
    fasta_path = Path(config["fasta"])
    gff_path = Path(config["gff"])
    snp_path = Path(config["snp_sites"])
    record = read_single_fasta(fasta_path)
    genome_len = len(record.seq)

    gff = Gff(gff_path, target_seqid=str(config["seqid"]))
    plot_len = max(genome_len, gff.genome_length)
    full_range = (0, plot_len)
    cds_plus = gff.extract_features("CDS", target_strand=1, target_range=full_range)
    cds_minus = gff.extract_features("CDS", target_strand=-1, target_range=full_range)
    trna = gff.extract_features("tRNA", target_range=full_range)
    rrna = gff.extract_features("rRNA", target_range=full_range)
    genes = gff.extract_features("gene", target_range=full_range)

    circos = Circos({label: plot_len}, start=0, end=360)
    sector = circos.sectors[0]

    outer = sector.add_track((96, 100))
    outer.axis(fc="#f7f7f7", ec="#333333", lw=0.6)
    outer.xticks_by_interval(
        int(config["tick_interval"]),
        label_formatter=lambda v: f"{int(v / 1000)} kb",
        label_size=7,
        line_kws={"lw": 0.6},
    )

    plus_track = sector.add_track((82, 94))
    plus_track.axis(fc="#ffffff", ec="#dddddd", lw=0.4)
    plus_track.genomic_features(
        cds_plus,
        plotstyle="arrow",
        r_lim=(83, 93),
        facecolor_handler=color_for_cds,
        lw=0.25,
        ec="#222222",
    )

    minus_track = sector.add_track((68, 80))
    minus_track.axis(fc="#ffffff", ec="#dddddd", lw=0.4)
    minus_track.genomic_features(
        cds_minus,
        plotstyle="arrow",
        r_lim=(69, 79),
        facecolor_handler=color_for_cds,
        lw=0.25,
        ec="#222222",
    )

    rna_track = sector.add_track((55, 66))
    rna_track.axis(fc="#ffffff", ec="#dddddd", lw=0.4)
    rna_track.genomic_features(
        trna,
        plotstyle="box",
        r_lim=(57, 61),
        facecolor_handler=lambda _: "#16a085",
        lw=0.2,
        ec="#222222",
    )
    rna_track.genomic_features(
        rrna,
        plotstyle="box",
        r_lim=(61.5, 65),
        facecolor_handler=lambda _: "#8e44ad",
        lw=0.2,
        ec="#222222",
    )

    gc_track = sector.add_track((39, 52))
    gc_track.axis(fc="#ffffff", ec="#dddddd", lw=0.4)
    window = 2_000 if genome_len < 180_000 else 5_000
    positions, gc_delta, genome_gc = sliding_gc(str(record.seq), window)
    gc_track.fill_between(
        positions,
        gc_delta,
        0,
        vmin=-0.25,
        vmax=0.25,
        color="#2f80ed",
        alpha=0.35,
    )
    gc_track.line(positions, gc_delta, vmin=-0.25, vmax=0.25, color="#2f80ed", lw=0.7)

    snp_positions = read_snp_positions(snp_path)
    snp_bin = 2_000 if label == "cpDNA" else 5_000
    snp_x, snp_counts = snp_density(snp_positions, plot_len, snp_bin)
    snp_track = sector.add_track((24, 36))
    snp_track.axis(fc="#ffffff", ec="#dddddd", lw=0.4)
    snp_track.bar(
        snp_x,
        snp_counts,
        width=snp_bin * 0.92,
        vmin=0,
        vmax=max(snp_counts) if snp_counts else 1,
        color="#c0392b",
        alpha=0.75,
        lw=0,
    )

    circos.text(str(config["title"]), r=20, size=12, weight="bold")
    if plot_len == genome_len:
        length_label = f"{genome_len:,} bp | GC {genome_gc:.1%}"
    else:
        length_label = f"FASTA {genome_len:,} bp | annotation span {plot_len:,} bp | GC {genome_gc:.1%}"
    circos.text(length_label, r=12, size=8, color="#444444")
    circos.text("+ strand CDS", r=88, deg=275, size=7, color="#444444")
    circos.text("- strand CDS", r=74, deg=275, size=7, color="#444444")
    circos.text("tRNA/rRNA", r=61, deg=275, size=7, color="#444444")
    circos.text("GC delta", r=45, deg=275, size=7, color="#444444")
    circos.text("SNP density", r=30, deg=275, size=7, color="#444444")

    fig = circos.plotfig(figsize=(11, 11), dpi=220)
    fig.patch.set_facecolor("white")
    legend_handles = [
        Patch(facecolor="#2f80ed", label="Photosynthesis / plastid genes"),
        Patch(facecolor="#27ae60", label="ATP synthase / NADH genes"),
        Patch(facecolor="#9b59b6", label="Ribosomal genes"),
        Patch(facecolor="#f2994a", label="RNA polymerase / matK"),
        Patch(facecolor="#d35400", label="Mitochondrial respiration genes"),
        Patch(facecolor="#16a085", label="tRNA"),
        Patch(facecolor="#8e44ad", label="rRNA"),
        Patch(facecolor="#b8b8b8", label="Lower-confidence calls"),
        Patch(facecolor="#c0392b", label="SNP density"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.01),
        ncol=4,
        fontsize=7,
        frameon=False,
    )
    for ext in ["png", "svg", "pdf"]:
        out_path = OUT_DIR / f"{label}.circular_genome_map.{ext}"
        fig.savefig(out_path, bbox_inches="tight", dpi=220)
    plt.close(fig)

    return {
        "organelle": label,
        "fasta": str(fasta_path.relative_to(ROOT)),
        "gff": str(gff_path.relative_to(ROOT)),
        "snp_sites": str(snp_path.relative_to(ROOT)),
        "length_bp": genome_len,
        "plot_span_bp": plot_len,
        "gc_fraction": f"{genome_gc:.5f}",
        "gene_features": len(genes),
        "cds_features": len(cds_plus) + len(cds_minus),
        "trna_features": len(trna),
        "rrna_features": len(rrna),
        "snp_sites_count": len(snp_positions),
        "snp_density_bin_bp": snp_bin,
        "png": f"{label}.circular_genome_map.png",
        "svg": f"{label}.circular_genome_map.svg",
        "pdf": f"{label}.circular_genome_map.pdf",
    }


def write_summary(rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "organelle",
        "fasta",
        "gff",
        "snp_sites",
        "length_bp",
        "plot_span_bp",
        "gc_fraction",
        "gene_features",
        "cds_features",
        "trna_features",
        "rrna_features",
        "snp_sites_count",
        "snp_density_bin_bp",
        "png",
        "svg",
        "pdf",
    ]
    with (OUT_DIR / "genome_map_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = [build_map(label, config) for label, config in INPUTS.items()]
    write_summary(rows)


if __name__ == "__main__":
    main()
