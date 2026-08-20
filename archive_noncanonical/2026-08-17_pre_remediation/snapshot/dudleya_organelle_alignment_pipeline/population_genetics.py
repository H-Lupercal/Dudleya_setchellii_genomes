"""Compute cpDNA and mtDNA population genetic summaries.

This stage consumes filtered haploid SNP alignments
and sample metadata, then writes pairwise Fst and per-population summary tables
for cpDNA and mtDNA separately.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

from dudleya_organelle_alignment_pipeline.pca_analysis import read_fasta, read_sample_metadata
from dudleya_organelle_alignment_pipeline.variant_calling import labeled_output_name


DEFAULT_SNP_ALIGNMENT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/10_snp_alignment"
)
DEFAULT_METADATA_PATH = Path(
    "dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/included_samples.tsv"
)
DEFAULT_OUTPUT_DIR = Path(
    "dudleya_organelle_alignment_pipeline/results/17_population_genetics"
)
DEFAULT_RUN_LABEL = "primary"
BASES = {"A", "C", "G", "T"}


class PopulationGeneticsError(RuntimeError):
    """Raised when population-genetic summaries cannot be computed."""


@dataclass(frozen=True)
class PopulationInput:
    organelle: str
    track_id: str
    sample_count: int
    alignment_sites: int
    missing_bases: int
    alignment_fasta_path: Path
    site_table_path: Path


@dataclass(frozen=True)
class PopulationResult:
    organelle: str
    track_id: str
    sample_count: int
    alignment_sites: int
    population_count: int
    pairwise_comparison_count: int
    pairwise_fst_path: Path
    population_summary_path: Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_population_inputs(
    snp_alignment_dir: Path = DEFAULT_SNP_ALIGNMENT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
) -> list[PopulationInput]:
    summary_path = snp_alignment_dir / labeled_output_name(
        "snp_alignment_summary.tsv",
        run_label,
    )
    inputs: list[PopulationInput] = []
    for row in read_tsv(summary_path):
        fasta_path = Path(row["alignment_fasta_path"])
        site_table_path = Path(row["site_table_path"])
        if not fasta_path.exists() or fasta_path.stat().st_size == 0:
            raise PopulationGeneticsError(f"Missing SNP alignment FASTA: {fasta_path}")
        inputs.append(
            PopulationInput(
                organelle=row["organelle"],
                track_id=row["track_id"],
                sample_count=int(row["sample_count"]),
                alignment_sites=int(row["alignment_sites"]),
                missing_bases=int(row["missing_bases"]),
                alignment_fasta_path=fasta_path,
                site_table_path=site_table_path,
            )
        )
    if not inputs:
        raise PopulationGeneticsError(f"No SNP alignment rows found in {summary_path}")
    return inputs


def population_code_for_sample(row: dict[str, str]) -> str:
    popcode = row.get("popcode", "").strip()
    if popcode:
        return popcode
    return ""


def group_sequences_by_population(
    records: list[tuple[str, str]],
    metadata: dict[str, dict[str, str]],
) -> dict[str, list[tuple[str, str]]]:
    groups: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for sample_id, sequence in records:
        popcode = population_code_for_sample(metadata.get(sample_id, {}))
        if popcode:
            groups[popcode].append((sample_id, sequence))
    return dict(sorted(groups.items()))


def allele_counts_at_site(sequences: list[str], site_index: int) -> Counter[str]:
    return Counter(
        sequence[site_index]
        for sequence in sequences
        if site_index < len(sequence) and sequence[site_index] in BASES
    )


def gene_diversity(counts: Counter[str]) -> float:
    total = sum(counts.values())
    if total <= 1:
        return 0.0
    return 1.0 - sum((count / total) ** 2 for count in counts.values())


def compute_pairwise_fst(pop1_sequences: list[str], pop2_sequences: list[str]) -> tuple[float, int]:
    if not pop1_sequences or not pop2_sequences:
        return 0.0, 0
    site_count = min(len(sequence) for sequence in [*pop1_sequences, *pop2_sequences])
    numerator = 0.0
    denominator = 0.0
    informative_sites = 0
    for site_index in range(site_count):
        counts1 = allele_counts_at_site(pop1_sequences, site_index)
        counts2 = allele_counts_at_site(pop2_sequences, site_index)
        if not counts1 or not counts2:
            continue
        combined = counts1 + counts2
        if len(combined) < 2:
            continue
        h1 = gene_diversity(counts1)
        h2 = gene_diversity(counts2)
        ht = gene_diversity(combined)
        if ht <= 0:
            continue
        hs = (h1 + h2) / 2.0
        numerator += max(ht - hs, 0.0)
        denominator += ht
        informative_sites += 1
    if denominator <= 0:
        return 0.0, informative_sites
    return numerator / denominator, informative_sites


def compute_haplotype_diversity(sequences: list[str]) -> float:
    n = len(sequences)
    if n <= 1:
        return 0.0
    counts = Counter(sequences)
    return (n / (n - 1)) * (1.0 - sum((count / n) ** 2 for count in counts.values()))


def compute_nucleotide_diversity(sequences: list[str]) -> float:
    if len(sequences) <= 1:
        return 0.0
    differences = 0
    compared_sites = 0
    for seq1, seq2 in combinations(sequences, 2):
        for base1, base2 in zip(seq1, seq2):
            if base1 not in BASES or base2 not in BASES:
                continue
            compared_sites += 1
            if base1 != base2:
                differences += 1
    if compared_sites == 0:
        return 0.0
    return differences / compared_sites


def private_variant_count(
    population: str,
    groups: dict[str, list[tuple[str, str]]],
    site_count: int,
) -> int:
    private_count = 0
    pop_sequences = [sequence for _, sequence in groups[population]]
    other_sequences = [
        sequence
        for other_population, records in groups.items()
        if other_population != population
        for _, sequence in records
    ]
    for site_index in range(site_count):
        pop_alleles = set(allele_counts_at_site(pop_sequences, site_index))
        other_alleles = set(allele_counts_at_site(other_sequences, site_index))
        if pop_alleles - other_alleles:
            private_count += 1
    return private_count


def run_one_population_summary(
    organelle: str,
    track_id: str,
    alignment_sites: int,
    alignment_fasta_path: Path,
    metadata_path: Path,
    output_dir: Path,
    run_label: str = DEFAULT_RUN_LABEL,
) -> PopulationResult:
    output_dir.mkdir(parents=True, exist_ok=True)
    records = read_fasta(alignment_fasta_path)
    metadata = read_sample_metadata(metadata_path)
    groups = group_sequences_by_population(records, metadata)
    if len(groups) < 2:
        raise PopulationGeneticsError(
            f"Need at least two metadata-resolved populations for {organelle}"
        )

    pairwise_rows: list[dict[str, str]] = []
    for pop1, pop2 in combinations(sorted(groups), 2):
        sequences1 = [sequence for _, sequence in groups[pop1]]
        sequences2 = [sequence for _, sequence in groups[pop2]]
        fst, informative_sites = compute_pairwise_fst(sequences1, sequences2)
        pairwise_rows.append(
            {
                "organelle": organelle,
                "population_1": pop1,
                "population_2": pop2,
                "n_population_1": str(len(sequences1)),
                "n_population_2": str(len(sequences2)),
                "informative_sites": str(informative_sites),
                "fst": f"{fst:.8f}",
            }
        )

    population_rows: list[dict[str, str]] = []
    for population, pop_records in groups.items():
        sample_ids = [sample_id for sample_id, _ in pop_records]
        sequences = [sequence for _, sequence in pop_records]
        unique_haplotypes = len(set(sequences))
        population_rows.append(
            {
                "organelle": organelle,
                "population": population,
                "sample_count": str(len(sequences)),
                "sample_ids": ",".join(sample_ids),
                "snp_sites": str(alignment_sites),
                "haplotype_count": str(unique_haplotypes),
                "haplotype_diversity": f"{compute_haplotype_diversity(sequences):.8f}",
                "nucleotide_diversity": f"{compute_nucleotide_diversity(sequences):.8f}",
                "private_variant_sites": str(
                    private_variant_count(population, groups, alignment_sites)
                ),
            }
        )

    prefix = output_dir / f"{organelle}.{run_label}.population_genetics"
    pairwise_path = Path(f"{prefix}.pairwise_fst.tsv")
    population_path = Path(f"{prefix}.population_summary.tsv")
    write_tsv(
        pairwise_path,
        pairwise_rows,
        [
            "organelle",
            "population_1",
            "population_2",
            "n_population_1",
            "n_population_2",
            "informative_sites",
            "fst",
        ],
    )
    write_tsv(
        population_path,
        population_rows,
        [
            "organelle",
            "population",
            "sample_count",
            "sample_ids",
            "snp_sites",
            "haplotype_count",
            "haplotype_diversity",
            "nucleotide_diversity",
            "private_variant_sites",
        ],
    )
    return PopulationResult(
        organelle=organelle,
        track_id=track_id,
        sample_count=len(records),
        alignment_sites=alignment_sites,
        population_count=len(groups),
        pairwise_comparison_count=len(pairwise_rows),
        pairwise_fst_path=pairwise_path,
        population_summary_path=population_path,
    )


def run_population_genetics(
    snp_alignment_dir: Path = DEFAULT_SNP_ALIGNMENT_DIR,
    metadata_path: Path = DEFAULT_METADATA_PATH,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    run_label: str = DEFAULT_RUN_LABEL,
) -> list[PopulationResult]:
    inputs = read_population_inputs(snp_alignment_dir=snp_alignment_dir, run_label=run_label)
    results = [
        run_one_population_summary(
            organelle=population_input.organelle,
            track_id=population_input.track_id,
            alignment_sites=population_input.alignment_sites,
            alignment_fasta_path=population_input.alignment_fasta_path,
            metadata_path=metadata_path,
            output_dir=output_dir,
            run_label=run_label,
        )
        for population_input in inputs
    ]
    write_population_genetics_outputs(output_dir, results, run_label=run_label)
    return results


def write_population_genetics_outputs(
    output_dir: Path,
    results: list[PopulationResult],
    run_label: str,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        output_dir / labeled_output_name("population_genetics_summary.tsv", run_label),
        [
            {
                "organelle": result.organelle,
                "track_id": result.track_id,
                "sample_count": str(result.sample_count),
                "alignment_sites": str(result.alignment_sites),
                "population_count": str(result.population_count),
                "pairwise_comparison_count": str(result.pairwise_comparison_count),
                "pairwise_fst_path": result.pairwise_fst_path.as_posix(),
                "population_summary_path": result.population_summary_path.as_posix(),
            }
            for result in results
        ],
        [
            "organelle",
            "track_id",
            "sample_count",
            "alignment_sites",
            "population_count",
            "pairwise_comparison_count",
            "pairwise_fst_path",
            "population_summary_path",
        ],
    )
    write_population_genetics_report(
        output_dir / labeled_output_name("population_genetics_report.md", run_label),
        results=results,
        run_label=run_label,
    )


def write_population_genetics_report(
    path: Path,
    results: list[PopulationResult],
    run_label: str,
) -> None:
    label = run_label or "full"
    lines = [
        "# Population Genetics",
        "",
        "This step computes pairwise population Fst and per-population summary",
        "statistics from filtered haploid cpDNA and mtDNA SNP alignments.",
        "Only samples with resolved population codes are included in these",
        "population-level summaries.",
        "",
        "## Run",
        "",
        f"- Run label: `{label}`",
        "- Fst: Nei-style haploid SNP differentiation, averaged across informative sites",
        "- Population summaries: sample count, haplotypes, diversity, nucleotide diversity, private variant sites",
        "",
        "## Results",
        "",
    ]
    for result in results:
        lines.extend(
            [
                f"### {result.organelle}",
                "",
                f"- Track: `{result.track_id}`",
                f"- Total samples in SNP alignment: {result.sample_count}",
                f"- Metadata-resolved populations: {result.population_count}",
                f"- Pairwise comparisons: {result.pairwise_comparison_count}",
                f"- Pairwise Fst table: `{result.pairwise_fst_path}`",
                f"- Population summary table: `{result.population_summary_path}`",
                "",
            ]
        )
    path.write_text("\n".join(lines))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute cpDNA/mtDNA Fst and population summaries."
    )
    parser.add_argument("--snp-alignment-dir", type=Path, default=DEFAULT_SNP_ALIGNMENT_DIR)
    parser.add_argument("--metadata-path", type=Path, default=DEFAULT_METADATA_PATH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-label", default=DEFAULT_RUN_LABEL)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    results = run_population_genetics(
        snp_alignment_dir=args.snp_alignment_dir,
        metadata_path=args.metadata_path,
        output_dir=args.output_dir,
        run_label=args.run_label,
    )
    for result in results:
        print(
            f"{result.organelle}: {result.population_count} populations, "
            f"{result.pairwise_comparison_count} pairwise Fst comparisons at "
            f"{result.pairwise_fst_path}"
        )
    print(f"Output directory: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
