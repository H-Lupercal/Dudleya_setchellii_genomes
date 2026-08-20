import os
import subprocess
from pathlib import Path

import pytest
from organelle_pipeline.variants import (
    Genotype,
    build_all_sites_call_command,
    build_primary_filter_commands,
    consensus_base,
    site_passes_filters,
)


def test_low_depth_or_quality_is_masked_not_replaced_by_reference() -> None:
    assert consensus_base("A", Genotype("T", depth=4, quality=99)) == "N"
    assert consensus_base("A", Genotype("T", depth=10, quality=19)) == "N"
    assert consensus_base("A", Genotype("T", depth=10, quality=20)) == "T"
    assert consensus_base("A", Genotype("A", depth=10, quality=20)) == "A"


def test_singletons_are_retained_for_primary_but_not_mac2_ordination() -> None:
    genotypes = ("A", "A", "A", "T")
    assert site_passes_filters(genotypes, quality=30, min_mac=1)
    assert not site_passes_filters(genotypes, quality=30, min_mac=2)


def test_all_site_call_and_filters_mask_genotypes_before_missingness() -> None:
    call = build_all_sites_call_command(
        "ref.fa",
        "samples.bam.list",
        "sites.bed",
        "all_sites.bcf",
        likelihood_bcf="likelihoods.bcf",
    )
    filters = build_primary_filter_commands(
        "ref.fa",
        "all_sites.bcf",
        "masked.bcf",
        "high_confidence.vcf.gz",
        "primary.vcf.gz",
        "mac2.vcf.gz",
    )
    rendered = "\n".join((call, *filters))
    assert "bcftools call --ploidy 1 -m -a FORMAT/GQ" in call
    assert "-Ob -o likelihoods.bcf" in call
    assert "FORMAT/DP,FORMAT/AD" in call
    assert "likelihoods.bcf &&" in call
    assert " -v " not in call
    assert "-d 250" in call
    assert "-q 20" in call
    assert "-Q 20" in call
    assert "Type=Integer" in call
    assert "Type=Float" in call
    assert "bcftools reheader" in call
    assert "+setGT" in filters[0]
    assert 'FMT/DP="." | FMT/GQ="." | FMT/DP<5 | FMT/GQ<20' in filters[0]
    assert "F_MISSING<=0.2" in rendered
    assert "MAC>=1" not in filters[1]
    assert "high_confidence.vcf.gz" in filters[1]
    assert "AC>=1 && AC<=AN-1" in rendered
    assert "AC>=2 && AC<=AN-2" in rendered


def test_bcftools_masks_only_the_failing_genotype(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[3]
    tool_dir = root / ".tools/bioconda-env/bin"
    if not (tool_dir / "bcftools").is_file():
        pytest.skip("pinned bcftools environment is unavailable")
    source = root / "canonical_publication/pipeline/tests/fixtures/genotype_filter_input.vcf"
    masked = tmp_path / "masked.bcf"
    command = build_primary_filter_commands(
        "unused.fa",
        source,
        masked,
        tmp_path / "unused.high_confidence.vcf.gz",
        tmp_path / "unused.primary.vcf.gz",
        tmp_path / "unused.mac2.vcf.gz",
    )[0]
    environment = os.environ.copy()
    environment["PATH"] = f"{tool_dir}{os.pathsep}{environment['PATH']}"

    subprocess.run(["bash", "-o", "pipefail", "-c", command], check=True, env=environment)
    observed = subprocess.run(
        [str(tool_dir / "bcftools"), "query", "-f", "%POS[\t%GT]\n", str(masked)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()

    assert observed == ["1\t1\t.", "2\t0\t.", "3\t.\t."]


def test_fixed_alternate_site_is_retained_for_consensus_but_not_segregating_sets(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[3]
    tool_dir = root / ".tools/bioconda-env/bin"
    if not all((tool_dir / executable).is_file() for executable in ("bcftools", "samtools")):
        pytest.skip("pinned bcftools/samtools environment is unavailable")
    reference = tmp_path / "reference.fa"
    source = tmp_path / "input.vcf"
    reference.write_text(">cp\nAC\n")
    source.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=cp,length=2>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\ts5\n"
        "cp\t1\t.\tA\tT\t40\tPASS\t.\tGT:DP:GQ\t1:10:30\t0:10:30\t0:10:30\t0:10:30\t0:10:30\n"
        "cp\t2\t.\tC\tG\t40\tPASS\t.\tGT:DP:GQ\t1:10:30\t1:10:30\t1:10:30\t1:10:30\t1:10:30\n"
    )
    masked = tmp_path / "masked.bcf"
    high_confidence = tmp_path / "high_confidence.vcf.gz"
    primary = tmp_path / "primary.vcf.gz"
    mac2 = tmp_path / "mac2.vcf.gz"
    commands = build_primary_filter_commands(
        reference,
        source,
        masked,
        high_confidence,
        primary,
        mac2,
    )
    environment = os.environ.copy()
    environment["PATH"] = f"{tool_dir}{os.pathsep}{environment['PATH']}"
    subprocess.run([str(tool_dir / "samtools"), "faidx", str(reference)], check=True)
    for command in commands:
        subprocess.run(
            ["bash", "-o", "pipefail", "-c", command],
            check=True,
            env=environment,
        )

    def positions(vcf: Path) -> list[str]:
        return subprocess.run(
            [str(tool_dir / "bcftools"), "query", "-f", "%POS\\n", str(vcf)],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()

    assert positions(high_confidence) == ["1", "2"]
    assert positions(primary) == ["1"]
    assert positions(mac2) == []
