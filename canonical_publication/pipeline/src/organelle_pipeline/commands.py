"""Shell command construction for external bioinformatics tools."""

from __future__ import annotations

import shlex
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class MappingInputs:
    sample_id: str
    r1_paths: tuple[str | Path, ...]
    r2_paths: tuple[str | Path, ...]
    reference: str | Path
    output_bam: str | Path
    fastp_json: str | Path
    fastp_html: str | Path


def _q(value: str | Path) -> str:
    return shlex.quote(str(value))


def build_mapping_command(
    inputs: MappingInputs,
    threads: int,
    qualified_quality_phred: int = 20,
    maximum_unqualified_base_percent: int = 40,
    minimum_length: int = 50,
    minimum_mapping_quality: int = 20,
    exclude_sam_flags: int = 3844,
    detect_adapters_for_pe: bool = True,
) -> str:
    """Stream fastp into BWA and produce a filtered, deduplicated BAM."""

    if threads < 2:
        raise ValueError("mapping requires at least two threads")
    if not inputs.r1_paths or len(inputs.r1_paths) != len(inputs.r2_paths):
        raise ValueError("mapping requires balanced paired-end lanes")
    if len(inputs.r1_paths) != 1:
        raise ValueError("streaming mapper currently requires one balanced lane per sample")
    output = Path(inputs.output_bam)
    sort_temp_dir = output.parent / "sort_tmp"
    name_sort_prefix = sort_temp_dir / f"{inputs.sample_id}.name"
    coordinate_sort_prefix = sort_temp_dir / f"{inputs.sample_id}.coordinate"
    rg = f"@RG\\tID:{inputs.sample_id}\\tSM:{inputs.sample_id}\\tPL:ILLUMINA"
    sort_threads = max(1, threads // 2)
    adapter_flag = " --detect_adapter_for_pe" if detect_adapters_for_pe else ""
    return (
        f"fastp --in1 {_q(inputs.r1_paths[0])} --in2 {_q(inputs.r2_paths[0])} "
        f"--qualified_quality_phred {qualified_quality_phred} "
        f"--unqualified_percent_limit {maximum_unqualified_base_percent} "
        f"--length_required {minimum_length}"
        f"{adapter_flag} --stdout "
        f"--json {_q(inputs.fastp_json)} --html {_q(inputs.fastp_html)} | "
        f"bwa mem -p -t {threads} -R {_q(rg)} {_q(inputs.reference)} - | "
        f"samtools view -F {exclude_sam_flags} -q {minimum_mapping_quality} -u - | "
        f"samtools sort -T {_q(name_sort_prefix)} -n -@ {sort_threads} -u - | "
        "samtools fixmate -m - - | "
        f"samtools sort -T {_q(coordinate_sort_prefix)} -@ {sort_threads} -u - | "
        f"samtools markdup -r -@ {sort_threads} - {_q(output)} && "
        f"samtools index {_q(output)}"
    )
