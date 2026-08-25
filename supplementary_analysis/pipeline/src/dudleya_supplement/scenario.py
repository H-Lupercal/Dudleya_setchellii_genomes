"""Threshold-sensitivity analyses sourced from immutable canonical mappings."""

from __future__ import annotations

import math
import os
import shutil
import subprocess
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
from Bio import SeqIO
from organelle_pipeline.analysis import build_iqtree_command
from organelle_pipeline.consensus import build_callable_sequence, reference_concordance
from organelle_pipeline.haplotypes import summarize_haplotypes
from organelle_pipeline.ordination import prepare_haploid_pca_matrix
from organelle_pipeline.popgen import callable_nucleotide_diversity, haplotype_diversity_from_assignments, hudson_fst
from organelle_pipeline.variants import build_all_sites_call_command, build_primary_filter_commands
from sklearn.decomposition import PCA

from .io import read_tsv, write_tsv
from .sensitivity import classify_fst, classify_pi, compare_pi, procrustes_permutation_test, rank_extreme_cases

ORGANELLES = ("chloroplast", "mitochondria")


def _fasta_records(path: Path) -> dict[str, str]:
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def _mask_array(path: Path, record: str, length: int, included: bool = True) -> np.ndarray:
    mask = np.zeros(length, dtype=bool) if included else np.ones(length, dtype=bool)
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        observed, start, end, *_ = line.split("\t")
        if observed == record:
            mask[int(start) : int(end)] = included
    return mask


def eligible_samples(root: Path, scenario: dict[str, object], organelle: str) -> list[str]:
    breadth = read_tsv(root / "canonical_publication/results/qc/publication-20260817/sample_breadth.tsv")
    prefix = "cp" if organelle == "chloroplast" else "mt"
    key = f"{prefix}_unique_sites_breadth_dp{int(scenario['eligibility_dp'])}"
    return sorted(row["sample_id"] for row in breadth if float(row[key]) >= float(scenario["breadth"]))


def write_scenario_metadata(root: Path, run_id: str, scenario_name: str, samples: dict[str, list[str]]) -> Path:
    corrected = read_tsv(root / "supplementary_analysis/metadata/samples/samples.corrected-20260824.tsv")
    by_sample = {row["sample_id"]: row for row in corrected}
    output = root / f"supplementary_analysis/results/sensitivity/{run_id}/{scenario_name}/metadata"
    output.mkdir(parents=True, exist_ok=True)
    for organelle in ORGANELLES:
        rows = [by_sample[sample] for sample in samples[organelle]]
        write_tsv(output / f"{organelle}_samples.tsv", rows, list(rows[0]), root)
    shared = sorted(set(samples["chloroplast"]) & set(samples["mitochondria"]))
    rows = [by_sample[sample] for sample in shared]
    write_tsv(output / "shared_samples.tsv", rows, list(rows[0]), root)
    return output


def _depth_cache(root: Path, run_id: str, sample: str) -> np.ndarray:
    cache = root / f"supplementary_analysis/work/{run_id}/depth/mitochondria/{sample}.npy"
    if cache.is_file():
        return np.load(cache)
    cache.parent.mkdir(parents=True, exist_ok=True)
    bam = root / f"canonical_publication/work/publication-20260817/mapping/{sample}.organelle.bam"
    process = subprocess.Popen(
        ["samtools", "depth", "-aa", "-q", "20", "-Q", "20", "-r", "mitochondria", str(bam)],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    depths = np.zeros(243359, dtype=np.uint16)
    for line in process.stdout:
        _, position, depth = line.rstrip().split("\t")
        depths[int(position) - 1] = min(int(depth), np.iinfo(np.uint16).max)
    if process.wait() != 0:
        raise RuntimeError(f"samtools depth failed for {sample}")
    np.save(cache, depths)
    return depths


def build_mitochondria_mask(
    root: Path,
    run_id: str,
    scenario_name: str,
    samples: list[str],
    *,
    depth: int,
    support_fraction: float,
) -> Path:
    repeat_mask = _mask_array(
        root / "canonical_publication/references/masks/mitochondria_repeat_mask.bed", "mitochondria", 243359, included=False
    )
    support = np.zeros(243359, dtype=np.uint16)
    for sample in samples:
        support += (_depth_cache(root, run_id, sample) >= depth).astype(np.uint16)
    accepted = repeat_mask & (support >= math.ceil(len(samples) * support_fraction))
    output = root / f"supplementary_analysis/results/sensitivity/{run_id}/{scenario_name}/masks/mitochondria_high_confidence_sites.bed"
    output.parent.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[int, int]] = []
    start: int | None = None
    for index, value in enumerate(np.append(accepted, False)):
        if value and start is None:
            start = index
        elif not value and start is not None:
            if index - start >= 50:
                rows.append((start, index))
            start = None
    with output.open("w") as handle:
        for start, end in rows:
            handle.write(f"mitochondria\t{start}\t{end}\tsupport_{support_fraction:.2f}\n")
    if not rows:
        raise RuntimeError(f"No mitochondrial mask intervals for {scenario_name}")
    return output


def call_scenario_variants(
    root: Path,
    run_id: str,
    scenario_name: str,
    organelle: str,
    samples: list[str],
    mask: Path,
    *,
    minimum_depth: int,
    minimum_gq: int,
    maximum_missing: float,
) -> dict[str, Path]:
    work = root / f"supplementary_analysis/work/{run_id}/sensitivity/{scenario_name}/variants"
    results = root / f"supplementary_analysis/results/sensitivity/{run_id}/{scenario_name}/variants"
    work.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    bam_list = work / f"{organelle}.bam.list"
    bam_list.write_text("".join(f"canonical_publication/work/publication-20260817/mapping/{sample}.organelle.bam\n" for sample in samples))
    likelihoods = work / f"{organelle}.mpileup_likelihoods.bcf"
    all_sites = work / f"{organelle}.all_sites.bcf"
    masked = work / f"{organelle}.genotype_masked.all_sites.bcf"
    high = results / f"{organelle}.high_confidence_variant_sites.vcf.gz"
    primary = results / f"{organelle}.primary.vcf.gz"
    mac2 = results / f"{organelle}.mac2_ordination.vcf.gz"
    if not mac2.is_file():
        commands = (
            build_all_sites_call_command(
                "canonical_publication/references/selected/organelle_combined.fa",
                bam_list.relative_to(root),
                mask.relative_to(root),
                all_sites.relative_to(root),
                likelihood_bcf=likelihoods.relative_to(root),
                minimum_mapping_quality=20,
                minimum_base_quality=20,
                maximum_per_file_depth=250,
                ploidy=1,
            ),
            *build_primary_filter_commands(
                "canonical_publication/references/selected/organelle_combined.fa",
                all_sites.relative_to(root),
                masked.relative_to(root),
                high.relative_to(root),
                primary.relative_to(root),
                mac2.relative_to(root),
                minimum_depth=minimum_depth,
                minimum_genotype_quality=minimum_gq,
                minimum_site_quality=30,
                maximum_missing_fraction=maximum_missing,
                primary_minimum_mac=1,
                ordination_minimum_mac=2,
            ),
        )
        log = work / f"{organelle}.commands.log"
        with log.open("w") as handle:
            for command in commands:
                handle.write(f"COMMAND\t{command}\n")
                handle.flush()
                subprocess.run(
                    ["bash", "-o", "pipefail", "-c", command],
                    cwd=root,
                    env=os.environ.copy(),
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    check=True,
                )
    return {"likelihoods": likelihoods, "all_sites": all_sites, "masked": masked, "high": high, "primary": primary, "mac2": mac2}


def _accepted_positions(vcf: Path) -> set[int]:
    result = subprocess.run(["bcftools", "query", "-f", "%POS\n", str(vcf)], capture_output=True, text=True, check=True)
    return {int(value) for value in result.stdout.splitlines() if value}


def _query_sample_rows(path: Path, sample: str, likelihoods: bool) -> list[tuple[object, ...]]:
    fmt = "%POS\t%REF\t%ALT[\t%DP\t%PL]\n" if likelihoods else "%POS\t%REF\t%ALT[\t%GT]\n"
    process = subprocess.Popen(["bcftools", "query", "-s", sample, "-f", fmt, str(path)], stdout=subprocess.PIPE, text=True)
    assert process.stdout is not None
    rows: list[tuple[object, ...]] = []
    for line in process.stdout:
        fields = line.rstrip().split("\t")
        if likelihoods:
            position, ref, alt, depth, pl = fields
            rows.append(
                (
                    int(position),
                    ref,
                    alt,
                    int(depth) if depth not in {"", "."} else -1,
                    tuple(int(value) for value in pl.split(",")) if pl not in {"", "."} else (),
                )
            )
        else:
            position, ref, alt, genotype = fields
            rows.append((int(position), ref, alt, genotype))
    if process.wait() != 0:
        raise RuntimeError(f"bcftools query failed for {sample}: {path}")
    return rows


def build_scenario_alignment(
    root: Path,
    run_id: str,
    scenario_name: str,
    organelle: str,
    samples: list[str],
    variants: dict[str, Path],
    *,
    minimum_depth: int,
    minimum_gq: int,
) -> tuple[Path, Path]:
    reference = _fasta_records(root / f"canonical_publication/references/selected/{organelle}.fa")[organelle]
    accepted = _accepted_positions(variants["high"])
    records: dict[str, str] = {}
    summaries = []
    for sample in samples:
        sequence = build_callable_sequence(
            reference,
            _query_sample_rows(variants["masked"], sample, False),  # type: ignore[arg-type]
            accepted,
            invariant_rows=_query_sample_rows(variants["likelihoods"], sample, True),  # type: ignore[arg-type]
            minimum_depth=minimum_depth,
            minimum_genotype_quality=minimum_gq,
        )
        concordance = reference_concordance(reference, sequence)
        records[sample] = sequence
        summaries.append(
            {
                "sample_id": sample,
                "reference_length": len(sequence),
                "callable_bases": concordance.callable_bases,
                "callable_fraction": f"{concordance.callable_bases / len(sequence):.12g}",
                "reference_matches": concordance.reference_matches,
                "nonreference_bases": concordance.nonreference_bases,
                "callable_reference_identity": f"{concordance.identity:.12g}",
            }
        )
    output = root / f"supplementary_analysis/results/sensitivity/{run_id}/{scenario_name}/alignments/{organelle}.callable_alignment.fa"
    summary = output.with_suffix(".summary.tsv")
    _write_fasta(output, records)
    write_tsv(summary, summaries, list(summaries[0]), root)
    return output, summary


def run_scenario_pca(root: Path, run_id: str, scenario_name: str, organelle: str, vcf: Path, metadata: Path) -> tuple[Path, Path]:
    samples = subprocess.run(["bcftools", "query", "-l", str(vcf)], capture_output=True, text=True, check=True).stdout.splitlines()
    queried = subprocess.run(["bcftools", "query", "-f", "%POS[\t%GT]\n", str(vcf)], capture_output=True, text=True, check=True)
    markers = [[float(value) if value in {"0", "1"} else np.nan for value in line.split("\t")[1:]] for line in queried.stdout.splitlines()]
    if not markers:
        raise RuntimeError(f"No MAC>=2 markers for PCA: {scenario_name} {organelle}")
    matrix = prepare_haploid_pca_matrix(np.asarray(markers, dtype=float).T)
    count = min(10, matrix.shape[0] - 1, matrix.shape[1])
    model = PCA(n_components=count, svd_solver="full")
    coordinates = model.fit_transform(matrix)
    metadata_rows = {row["sample_id"]: row for row in read_tsv(metadata)}
    output_dir = root / f"supplementary_analysis/results/sensitivity/{run_id}/{scenario_name}/pca"
    output_dir.mkdir(parents=True, exist_ok=True)
    coordinate_path = output_dir / f"{organelle}.coordinates.tsv"
    variance_path = output_dir / f"{organelle}.variance.tsv"
    fields = ["sample_id", "popcode", *[f"PC{i}" for i in range(1, count + 1)]]
    rows = [
        {
            "sample_id": sample,
            "popcode": metadata_rows[sample]["popcode"],
            **{f"PC{i + 1}": f"{value:.12g}" for i, value in enumerate(values)},
        }
        for sample, values in zip(samples, coordinates, strict=True)
    ]
    write_tsv(coordinate_path, rows, fields, root)
    write_tsv(
        variance_path,
        [
            {"component": f"PC{i}", "explained_variance_ratio": f"{value:.12g}"}
            for i, value in enumerate(model.explained_variance_ratio_, 1)
        ],
        ["component", "explained_variance_ratio"],
        root,
    )
    return coordinate_path, variance_path


def run_scenario_haplotypes(
    root: Path, run_id: str, scenario_name: str, organelle: str, alignment: Path, metadata: Path
) -> dict[str, Path]:
    records = _fasta_records(alignment)
    summary = summarize_haplotypes(records)
    meta = {row["sample_id"]: row for row in read_tsv(metadata)}
    output_dir = root / f"supplementary_analysis/results/sensitivity/{run_id}/{scenario_name}/haplotypes"
    assignments = output_dir / f"{organelle}.sample_haplotypes.tsv"
    haplotypes = output_dir / f"{organelle}.haplotypes.tsv"
    positions = output_dir / f"{organelle}.haplotype_positions.tsv"
    edges = output_dir / f"{organelle}.network_edges.tsv"
    assignment_rows = [
        {"sample_id": sample, "popcode": meta[sample]["popcode"], "haplotype": summary.sample_haplotypes[sample]}
        for sample in sorted(records)
    ]
    write_tsv(assignments, assignment_rows, ["sample_id", "popcode", "haplotype"], root)
    hap_rows = []
    for haplotype in sorted(summary.counts):
        counts = Counter(meta[sample]["popcode"] for sample, assigned in summary.sample_haplotypes.items() if assigned == haplotype)
        hap_rows.append(
            {
                "haplotype": haplotype,
                "sample_count": summary.counts[haplotype],
                "population_counts": ",".join(f"{key}:{value}" for key, value in sorted(counts.items())),
                "sequence": summary.sequences[haplotype],
            }
        )
    write_tsv(haplotypes, hap_rows, ["haplotype", "sample_count", "population_counts", "sequence"], root)
    write_tsv(positions, [{"alignment_position_1based": value + 1} for value in summary.positions], ["alignment_position_1based"], root)
    graph = nx.Graph()
    graph.add_nodes_from(summary.counts)
    for left, right in combinations(sorted(summary.counts), 2):
        graph.add_edge(left, right, weight=sum(a != b for a, b in zip(summary.sequences[left], summary.sequences[right], strict=True)))
    tree = nx.minimum_spanning_tree(graph, weight="weight", algorithm="kruskal")
    write_tsv(
        edges,
        [
            {"haplotype_1": left, "haplotype_2": right, "mutational_distance": values["weight"]}
            for left, right, values in sorted(tree.edges(data=True))
        ],
        ["haplotype_1", "haplotype_2", "mutational_distance"],
        root,
    )
    return {"assignments": assignments, "haplotypes": haplotypes, "positions": positions, "edges": edges}


def run_scenario_popgen(
    root: Path,
    run_id: str,
    scenario_name: str,
    organelle: str,
    alignment: Path,
    metadata: Path,
    assignments: Path,
) -> tuple[Path, Path]:
    records = _fasta_records(alignment)
    meta = {row["sample_id"]: row for row in read_tsv(metadata)}
    hap = {row["sample_id"]: row["haplotype"] for row in read_tsv(assignments)}
    groups: dict[str, list[str]] = defaultdict(list)
    group_samples: dict[str, list[str]] = defaultdict(list)
    for sample, sequence in records.items():
        if meta[sample]["population_inference_eligible"] != "yes":
            continue
        groups[meta[sample]["popcode"]].append(sequence)
        group_samples[meta[sample]["popcode"]].append(sample)
    population_rows = []
    for population in sorted(groups):
        estimate = callable_nucleotide_diversity(groups[population])
        hap_estimate = haplotype_diversity_from_assignments([hap[sample] for sample in group_samples[population]])
        population_rows.append(
            {
                "organelle": organelle,
                "population": population,
                "sample_count": len(groups[population]),
                "sample_ids": ",".join(group_samples[population]),
                "pairwise_differences": estimate.differences,
                "jointly_callable_sites": estimate.jointly_callable_sites,
                "pairwise_callable_sites": estimate.compared_sites,
                "nucleotide_diversity": f"{estimate.pi:.12g}",
                "haplotype_count": hap_estimate.haplotype_count,
                "haplotype_diversity": f"{hap_estimate.diversity:.12g}",
                "haplotype_assigned_samples": hap_estimate.assigned_samples,
                "haplotype_ambiguous_samples": hap_estimate.ambiguous_samples,
            }
        )
    pair_rows = []
    for left, right in combinations(sorted(groups), 2):
        estimate = hudson_fst(groups[left], groups[right])
        pair_rows.append(
            {
                "organelle": organelle,
                "population_1": left,
                "population_2": right,
                "n_population_1": len(groups[left]),
                "n_population_2": len(groups[right]),
                "numerator": f"{estimate.numerator:.12g}",
                "denominator": f"{estimate.denominator:.12g}",
                "callable_sites_with_at_least_two_calls_per_population": estimate.callable_sites,
                "hudson_fst": f"{estimate.fst:.12g}",
            }
        )
    expected = len(groups) * (len(groups) - 1) // 2
    if len(pair_rows) != expected:
        raise RuntimeError("Population pair count invariant failed")
    output_dir = root / f"supplementary_analysis/results/sensitivity/{run_id}/{scenario_name}/popgen"
    population_path = output_dir / f"{organelle}.population_summary.tsv"
    pair_path = output_dir / f"{organelle}.pairwise_hudson_fst.tsv"
    write_tsv(population_path, population_rows, list(population_rows[0]), root)
    write_tsv(pair_path, pair_rows, list(pair_rows[0]), root)
    return population_path, pair_path


def run_scenario_tree(root: Path, run_id: str, scenario_name: str, organelle: str, alignment: Path, seed: int) -> Path:
    output_dir = root / f"supplementary_analysis/results/sensitivity/{run_id}/{scenario_name}/trees"
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / f"{organelle}.primary"
    treefile = Path(f"{prefix}.treefile")
    if treefile.is_file():
        return treefile
    command = (
        build_iqtree_command(
            alignment.relative_to(root),
            prefix.relative_to(root),
            seed=seed,
            model="MFP",
            sh_alrt_replicates=1000,
            ultrafast_bootstrap_replicates=1000,
            bnni=True,
        )
        + " -nt 8"
    )
    log = root / f"supplementary_analysis/work/{run_id}/sensitivity/{scenario_name}/trees/{organelle}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as handle:
        subprocess.run(["bash", "-o", "pipefail", "-c", command], cwd=root, stdout=handle, stderr=subprocess.STDOUT, check=True)
    return treefile


def run_full_scenario(root: Path, run_id: str, scenario_name: str, scenario: dict[str, object], seeds: dict[str, object]) -> list[Path]:
    samples = {organelle: eligible_samples(root, scenario, organelle) for organelle in ORGANELLES}
    metadata_dir = write_scenario_metadata(root, run_id, scenario_name, samples)
    mt_mask = build_mitochondria_mask(
        root,
        run_id,
        scenario_name,
        samples["mitochondria"],
        depth=int(scenario["eligibility_dp"]),
        support_fraction=0.80,
    )
    outputs: list[Path] = [
        mt_mask,
        metadata_dir / "chloroplast_samples.tsv",
        metadata_dir / "mitochondria_samples.tsv",
        metadata_dir / "shared_samples.tsv",
    ]
    for organelle in ORGANELLES:
        mask = root / "canonical_publication/references/masks/chloroplast_population_sites.bed" if organelle == "chloroplast" else mt_mask
        variants = call_scenario_variants(
            root,
            run_id,
            scenario_name,
            organelle,
            samples[organelle],
            mask,
            minimum_depth=int(scenario["dp"]),
            minimum_gq=int(scenario["gq"]),
            maximum_missing=float(scenario["missing"]),
        )
        alignment, summary = build_scenario_alignment(
            root,
            run_id,
            scenario_name,
            organelle,
            samples[organelle],
            variants,
            minimum_depth=int(scenario["dp"]),
            minimum_gq=int(scenario["gq"]),
        )
        pca = run_scenario_pca(root, run_id, scenario_name, organelle, variants["mac2"], metadata_dir / f"{organelle}_samples.tsv")
        haplotypes = run_scenario_haplotypes(root, run_id, scenario_name, organelle, alignment, metadata_dir / f"{organelle}_samples.tsv")
        popgen = run_scenario_popgen(
            root,
            run_id,
            scenario_name,
            organelle,
            alignment,
            metadata_dir / f"{organelle}_samples.tsv",
            haplotypes["assignments"],
        )
        tree = run_scenario_tree(
            root,
            run_id,
            scenario_name,
            organelle,
            alignment,
            int(seeds["cp_tree"] if organelle == "chloroplast" else seeds["mt_tree"]),
        )
        outputs.extend([*variants.values(), alignment, summary, *pca, *haplotypes.values(), *popgen, tree])
    return outputs


def _copy_canonical_scenario(root: Path, run_id: str) -> list[Path]:
    """Materialize a byte-copy baseline while recalculating DUSE-dependent summaries."""
    destination = root / f"supplementary_analysis/results/sensitivity/{run_id}/canonical"
    corrected = read_tsv(root / "supplementary_analysis/metadata/samples/samples.corrected-20260824.tsv")
    by_sample = {row["sample_id"]: row for row in corrected}
    outputs: list[Path] = []
    metadata_dir = destination / "metadata"
    for organelle in ORGANELLES:
        source_metadata = read_tsv(root / f"canonical_publication/metadata/qc/publication-20260817/{organelle}_samples.tsv")
        rows = [by_sample[row["sample_id"]] for row in source_metadata]
        metadata = metadata_dir / f"{organelle}_samples.tsv"
        write_tsv(metadata, rows, list(rows[0]), root)
        outputs.append(metadata)
        for category, filenames in {
            "variants": [f"{organelle}.primary.vcf.gz", f"{organelle}.mac2_ordination.vcf.gz"],
            "alignments": [f"{organelle}.callable_alignment.fa", f"{organelle}.callable_summary.tsv"],
            "pca": [f"{organelle}.coordinates.tsv", f"{organelle}.variance.tsv"],
            "haplotypes": [f"{organelle}.sample_haplotypes.tsv", f"{organelle}.haplotypes.tsv"],
            "trees": [f"{organelle}.primary.treefile"],
        }.items():
            canonical_category = "trees" if category == "trees" else category
            for filename in filenames:
                source = root / f"canonical_publication/results/{canonical_category}/publication-20260817/{filename}"
                target = destination / category / filename
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)
                outputs.append(target)
        popgen = run_scenario_popgen(
            root,
            run_id,
            "canonical",
            organelle,
            destination / f"alignments/{organelle}.callable_alignment.fa",
            metadata,
            destination / f"haplotypes/{organelle}.sample_haplotypes.tsv",
        )
        outputs.extend(popgen)
    shared_ids = sorted(
        set(row["sample_id"] for row in read_tsv(metadata_dir / "chloroplast_samples.tsv"))
        & set(row["sample_id"] for row in read_tsv(metadata_dir / "mitochondria_samples.tsv"))
    )
    shared = [by_sample[sample] for sample in shared_ids]
    shared_path = metadata_dir / "shared_samples.tsv"
    write_tsv(shared_path, shared, list(shared[0]), root)
    outputs.append(shared_path)
    return outputs


def run_mt_mask_scenario(
    root: Path,
    run_id: str,
    scenario_name: str,
    support_fraction: float,
    canonical_scenario: dict[str, object],
    seeds: dict[str, object],
) -> list[Path]:
    samples = eligible_samples(root, canonical_scenario, "mitochondria")
    metadata_dir = write_scenario_metadata(
        root,
        run_id,
        scenario_name,
        {"chloroplast": eligible_samples(root, canonical_scenario, "chloroplast"), "mitochondria": samples},
    )
    mask = build_mitochondria_mask(root, run_id, scenario_name, samples, depth=5, support_fraction=support_fraction)
    variants = call_scenario_variants(
        root,
        run_id,
        scenario_name,
        "mitochondria",
        samples,
        mask,
        minimum_depth=5,
        minimum_gq=20,
        maximum_missing=0.20,
    )
    alignment, alignment_summary = build_scenario_alignment(
        root, run_id, scenario_name, "mitochondria", samples, variants, minimum_depth=5, minimum_gq=20
    )
    metadata = metadata_dir / "mitochondria_samples.tsv"
    pca = run_scenario_pca(root, run_id, scenario_name, "mitochondria", variants["mac2"], metadata)
    haplotypes = run_scenario_haplotypes(root, run_id, scenario_name, "mitochondria", alignment, metadata)
    popgen = run_scenario_popgen(root, run_id, scenario_name, "mitochondria", alignment, metadata, haplotypes["assignments"])
    tree = run_scenario_tree(root, run_id, scenario_name, "mitochondria", alignment, int(seeds["mt_tree"]))
    return [
        mask,
        *(metadata_dir.glob("*.tsv")),
        *variants.values(),
        alignment,
        alignment_summary,
        *pca,
        *haplotypes.values(),
        *popgen,
        tree,
    ]


def _numeric_map(path: Path, key_fields: tuple[str, ...], value_field: str) -> dict[str, float]:
    result: dict[str, float] = {}
    for row in read_tsv(path):
        value = row[value_field]
        if value.lower() == "nan":
            continue
        result["|".join(row[field] for field in key_fields)] = float(value)
    return result


def _coordinate_matrix(path: Path, common: list[str]) -> np.ndarray:
    rows = {row["sample_id"]: row for row in read_tsv(path)}
    return np.asarray([[float(rows[sample][f"PC{component}"]) for component in range(1, 4)] for sample in common])


def summarize_sensitivity(root: Path, run_id: str, protest_seeds: list[int]) -> tuple[Path, Path, Path]:
    base = root / f"supplementary_analysis/results/sensitivity/{run_id}"
    scenarios = ("permissive", "strict", "mtmask70", "mtmask90")
    summary_rows: list[dict[str, object]] = []
    status_rows: list[dict[str, object]] = []
    extreme_rows: list[dict[str, object]] = []
    seed_index = 0
    for scenario_name in scenarios:
        organelles = ORGANELLES if scenario_name in {"permissive", "strict"} else ("mitochondria",)
        for organelle in organelles:
            canonical_pi = _numeric_map(
                base / f"canonical/popgen/{organelle}.population_summary.tsv", ("population",), "nucleotide_diversity"
            )
            scenario_pi = _numeric_map(
                base / f"{scenario_name}/popgen/{organelle}.population_summary.tsv", ("population",), "nucleotide_diversity"
            )
            pi = compare_pi(canonical_pi, scenario_pi)
            canonical_fst = _numeric_map(
                base / f"canonical/popgen/{organelle}.pairwise_hudson_fst.tsv", ("population_1", "population_2"), "hudson_fst"
            )
            scenario_fst = _numeric_map(
                base / f"{scenario_name}/popgen/{organelle}.pairwise_hudson_fst.tsv", ("population_1", "population_2"), "hudson_fst"
            )
            common_pairs = sorted(set(canonical_fst) & set(scenario_fst))
            left = [canonical_fst[key] for key in common_pairs]
            right = [scenario_fst[key] for key in common_pairs]
            from scipy.stats import spearmanr

            rho = float(spearmanr(left, right).statistic) if len(common_pairs) > 1 else float("nan")
            deltas = [abs(a - b) for a, b in zip(left, right, strict=True)]
            extreme_rows.extend(
                rank_extreme_cases(
                    canonical_pi,
                    scenario_pi,
                    scenario=scenario_name,
                    organelle=organelle,
                    metric="pi",
                )
            )
            extreme_rows.extend(
                rank_extreme_cases(
                    canonical_fst,
                    scenario_fst,
                    scenario=scenario_name,
                    organelle=organelle,
                    metric="fst",
                )
            )
            canonical_coords = {row["sample_id"] for row in read_tsv(base / f"canonical/pca/{organelle}.coordinates.tsv")}
            scenario_coords = {row["sample_id"] for row in read_tsv(base / f"{scenario_name}/pca/{organelle}.coordinates.tsv")}
            common_samples = sorted(canonical_coords & scenario_coords)
            seed = protest_seeds[seed_index]
            seed_index += 1
            protest = procrustes_permutation_test(
                _coordinate_matrix(base / f"canonical/pca/{organelle}.coordinates.tsv", common_samples),
                _coordinate_matrix(base / f"{scenario_name}/pca/{organelle}.coordinates.tsv", common_samples),
                permutations=9999,
                seed=seed,
            )
            primary_vcf = base / f"{scenario_name}/variants/{organelle}.primary.vcf.gz"
            snps = int(
                subprocess.run(["bcftools", "view", "-H", str(primary_vcf)], capture_output=True, text=True, check=True).stdout.count("\n")
            )
            summary_rows.append(
                {
                    "scenario": scenario_name,
                    "organelle": organelle,
                    "eligible_samples": len(read_tsv(base / f"{scenario_name}/metadata/{organelle}_samples.tsv")),
                    "retained_snps": snps,
                    "pi_common_populations": len(set(canonical_pi) & set(scenario_pi)),
                    "pi_spearman_rho": f"{pi.spearman_rho:.12g}",
                    "pi_median_proportional_change": f"{pi.median_proportional_change:.12g}",
                    "pi_maximum_proportional_change": f"{pi.maximum_proportional_change:.12g}",
                    "pi_zero_to_nonzero": pi.zero_to_nonzero,
                    "pi_nonzero_to_zero": pi.nonzero_to_zero,
                    "fst_common_pairs": len(common_pairs),
                    "fst_spearman_rho": f"{rho:.12g}",
                    "fst_median_absolute_change": f"{np.median(deltas):.12g}",
                    "fst_maximum_absolute_change": f"{max(deltas):.12g}",
                    "pca_common_samples": len(common_samples),
                    "protest_r": f"{protest.correlation:.12g}",
                    "protest_p": f"{protest.p_value:.12g}",
                    "protest_permutations": protest.permutations,
                    "protest_seed": protest.seed,
                }
            )
            status_rows.extend(
                [
                    {
                        "scenario": scenario_name,
                        "organelle": organelle,
                        "metric": "pi",
                        "status": classify_pi(pi.spearman_rho, pi.median_proportional_change),
                    },
                    {
                        "scenario": scenario_name,
                        "organelle": organelle,
                        "metric": "fst",
                        "status": classify_fst(rho, float(np.median(deltas))),
                    },
                    {
                        "scenario": scenario_name,
                        "organelle": organelle,
                        "metric": "pca",
                        "status": "PASS"
                        if protest.correlation >= 0.90 and protest.p_value < 0.001
                        else "PASS_WITH_CAVEAT"
                        if protest.correlation >= 0.80
                        else "FAIL",
                    },
                ]
            )
    summary_path = base / "sensitivity_summary.tsv"
    status_path = base / "sensitivity_status.tsv"
    extremes_path = base / "sensitivity_extreme_cases.tsv"
    write_tsv(summary_path, summary_rows, list(summary_rows[0]), root)
    write_tsv(status_path, status_rows, ["scenario", "organelle", "metric", "status"], root)
    write_tsv(
        extremes_path,
        extreme_rows,
        [
            "scenario",
            "organelle",
            "metric",
            "rank",
            "population_1",
            "population_2",
            "canonical_value",
            "scenario_value",
            "signed_change",
            "absolute_change",
            "proportional_change",
            "transition_type",
        ],
        root,
    )
    return summary_path, status_path, extremes_path


def run_all_sensitivity(root: Path, run_id: str, config: dict[str, object]) -> list[Path]:
    outputs = _copy_canonical_scenario(root, run_id)
    scenarios = config["scenarios"]  # type: ignore[index]
    seeds = config["seeds"]  # type: ignore[index]
    for name in ("permissive", "strict"):
        outputs.extend(run_full_scenario(root, run_id, name, scenarios[name], seeds))  # type: ignore[index]
    outputs.extend(run_mt_mask_scenario(root, run_id, "mtmask70", 0.70, scenarios["canonical"], seeds))  # type: ignore[index]
    outputs.extend(run_mt_mask_scenario(root, run_id, "mtmask90", 0.90, scenarios["canonical"], seeds))  # type: ignore[index]
    outputs.extend(summarize_sensitivity(root, run_id, [int(value) for value in seeds["protest"]]))  # type: ignore[index]
    return outputs
