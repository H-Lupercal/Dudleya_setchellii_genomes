"""Haploid variant filtering and callable-consensus rules."""

from __future__ import annotations

import shlex
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

BASES = frozenset("ACGT")


@dataclass(frozen=True)
class Genotype:
    allele: str | None
    depth: int
    quality: float | None


def consensus_base(
    reference: str,
    genotype: Genotype,
    min_depth: int = 5,
    min_quality: float = 20,
) -> str:
    """Return a callable haploid base or N when evidence is insufficient."""

    allele = (genotype.allele or "").upper()
    if genotype.depth < min_depth or genotype.quality is None or genotype.quality < min_quality or allele not in BASES:
        return "N"
    return allele


def site_passes_filters(
    genotypes: tuple[str, ...] | list[str],
    quality: float,
    min_mac: int = 1,
    minimum_quality: float = 30,
    maximum_missing: float = 0.20,
) -> bool:
    """Apply biallelic haploid SNP, QUAL, missingness, and MAC filters."""

    if quality < minimum_quality or not genotypes:
        return False
    called = [value.upper() for value in genotypes if value.upper() in BASES]
    if 1 - (len(called) / len(genotypes)) > maximum_missing:
        return False
    counts = Counter(called)
    if len(counts) != 2:
        return False
    return min(counts.values()) >= min_mac


def _quote(value: str | Path) -> str:
    return shlex.quote(str(value))


def build_all_sites_call_command(
    reference: str | Path,
    bam_list: str | Path,
    regions: str | Path,
    output_bcf: str | Path,
    likelihood_bcf: str | Path | None = None,
    minimum_mapping_quality: int = 20,
    minimum_base_quality: int = 20,
    maximum_per_file_depth: int = 250,
    ploidy: int = 1,
) -> str:
    """Call haploid genotypes and repair bcftools 1.24's invalid MQ header.

    Bcftools 1.24 emits the reserved INFO/MQ field as ``Type=Integer`` even
    though the VCF specification requires ``Type=Float``.  Values emitted by
    the caller are integers and therefore remain valid floating-point values;
    only the declaration is repaired.  Keeping the repair in the recorded
    command prevents warning-bearing BCFs from reaching consensus generation.
    """

    likelihood_bcf = Path(likelihood_bcf) if likelihood_bcf is not None else Path(f"{output_bcf}.likelihoods.bcf")
    uncorrected_bcf = Path(f"{output_bcf}.uncorrected")
    corrected_header = Path(f"{output_bcf}.header")
    mq_integer = "##INFO=<ID=MQ,Number=1,Type=Integer"
    mq_float = "##INFO=<ID=MQ,Number=1,Type=Float"

    return (
        f"bcftools mpileup -f {_quote(reference)} -b {_quote(bam_list)} "
        f"-R {_quote(regions)} -d {maximum_per_file_depth} "
        f"-q {minimum_mapping_quality} -Q {minimum_base_quality} "
        f"-a FORMAT/DP,FORMAT/AD -Ob -o {_quote(likelihood_bcf)} && "
        f"bcftools call --ploidy {ploidy} -m -a FORMAT/GQ "
        f"-Ob -o {_quote(uncorrected_bcf)} {_quote(likelihood_bcf)} && "
        f"bcftools view -h {_quote(uncorrected_bcf)} | "
        f"sed 's/{mq_integer}/{mq_float}/' > {_quote(corrected_header)} && "
        f"bcftools reheader -h {_quote(corrected_header)} "
        f"-o {_quote(output_bcf)} {_quote(uncorrected_bcf)} && "
        f"bcftools view -h {_quote(output_bcf)} | grep -Fq '{mq_float}' && "
        f"rm -f {_quote(uncorrected_bcf)} {_quote(corrected_header)}"
    )


def build_primary_filter_commands(
    reference: str | Path,
    all_sites_bcf: str | Path,
    masked_bcf: str | Path,
    high_confidence_vcf: str | Path,
    primary_vcf: str | Path,
    mac2_vcf: str | Path,
    minimum_depth: int = 5,
    minimum_genotype_quality: int = 20,
    minimum_site_quality: int = 30,
    maximum_missing_fraction: float = 0.20,
    primary_minimum_mac: int = 1,
    ordination_minimum_mac: int = 2,
) -> tuple[str, ...]:
    """Mask low-confidence GTs, then apply site and analysis-specific filters."""

    mask = (
        f"bcftools +setGT {_quote(all_sites_bcf)} -Ob -o {_quote(masked_bcf)} -- "
        f'-t q -n . -i \'FMT/DP="." | FMT/GQ="." | '
        f"FMT/DP<{minimum_depth} | FMT/GQ<{minimum_genotype_quality}'"
    )
    high_confidence = (
        f"bcftools norm -f {_quote(reference)} -Ou {_quote(masked_bcf)} | "
        "bcftools view -m2 -M2 -v snps -Ou | "
        "bcftools +fill-tags -Ou -- -t AC,AN,F_MISSING | "
        f"bcftools view -i 'QUAL>={minimum_site_quality} && "
        f"F_MISSING<={maximum_missing_fraction:.12g}' "
        f"-Oz -o {_quote(high_confidence_vcf)}"
    )
    primary = (
        f"bcftools view -i 'AC>={primary_minimum_mac} && "
        f"AC<=AN-{primary_minimum_mac}' "
        f"-Oz -o {_quote(primary_vcf)} {_quote(high_confidence_vcf)}"
    )
    mac2 = (
        f"bcftools view -i 'AC>={ordination_minimum_mac} && "
        f"AC<=AN-{ordination_minimum_mac}' "
        f"-Oz -o {_quote(mac2_vcf)} {_quote(primary_vcf)}"
    )
    return (
        mask,
        high_confidence,
        f"bcftools index -f {_quote(high_confidence_vcf)}",
        primary,
        f"bcftools index -f {_quote(primary_vcf)}",
        mac2,
        f"bcftools index -f {_quote(mac2_vcf)}",
    )
