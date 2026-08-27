import ast
from pathlib import Path

from organelle_pipeline.commands import MappingInputs, build_mapping_command


def _pipeline_source() -> tuple[str, tuple[str, ...]]:
    root = Path(__file__).resolve().parents[3]
    source = (root / "canonical_publication/pipeline/scripts/run_pipeline.py").read_text()
    module = ast.parse(source)
    assignment = next(
        node
        for node in module.body
        if isinstance(node, ast.Assign) and any(getattr(target, "id", "") == "STAGES" for target in node.targets)
    )
    return source, ast.literal_eval(assignment.value)


def test_publication_figures_run_after_admixture_and_before_reports() -> None:
    source, stages = _pipeline_source()

    assert stages.index("admixture") < stages.index("figures") < stages.index("reports")
    assert '"figures": [' in source
    assert 'str(scripts / "render_figures.py")' in source
    assert '"--config",' in source


def test_pairwise_distances_run_after_consensus_and_before_figures_and_reports() -> None:
    source, stages = _pipeline_source()

    assert stages.index("consensus") < stages.index("distances") < stages.index("figures") < stages.index("reports")
    assert '"distances": [' in source
    assert 'str(scripts / "sample_distances.py")' in source


def test_report_acceptance_explicitly_requires_figure_state_and_manifest() -> None:
    root = Path(__file__).resolve().parents[3]
    source = (root / "canonical_publication/pipeline/scripts/build_reports.py").read_text()

    assert 'run_provenance_dir / "figures.json"' in source
    assert "figure_manifest.tsv" in source


def test_report_acceptance_requires_validated_sample_distance_outputs() -> None:
    root = Path(__file__).resolve().parents[3]
    source = (root / "canonical_publication/pipeline/scripts/build_reports.py").read_text()

    assert "validate_pairwise_distance_outputs" in source
    assert 'run_provenance_dir / "distances"' in source


def test_report_acceptance_distinguishes_provider_manifest_self_reference_from_read_failure() -> None:
    root = Path(__file__).resolve().parents[3]
    source = (root / "canonical_publication/pipeline/scripts/build_reports.py").read_text()

    assert '"UNVERIFIABLE_SELF_REFERENCE"' in source


def test_report_acceptance_scans_provenance_logs_for_absolute_paths() -> None:
    root = Path(__file__).resolve().parents[3]
    source = (root / "canonical_publication/pipeline/scripts/build_reports.py").read_text()

    assert 'if "/logs/" not in path' not in source


def test_report_fingerprint_tracks_the_full_manifested_pipeline_code_surface() -> None:
    root = Path(__file__).resolve().parents[3]
    source = (root / "canonical_publication/pipeline/scripts/build_reports.py").read_text()

    assert "pipeline_code_digest(root)" in source


def test_finalized_mapping_can_be_rebound_after_downstream_only_code_change() -> None:
    root = Path(__file__).resolve().parents[3]
    mapper = (root / "canonical_publication/pipeline/scripts/map_samples.py").read_text()
    finalizer = (root / "canonical_publication/pipeline/scripts/finalize_mapping_provenance.py").read_text()

    assert "pending provenance rebind" in mapper
    assert 'raise StaleOutputError(f"Finalized pipeline code is stale' not in mapper
    assert '"provenance_rebind"' in finalizer
    assert "rebound_from_fingerprint" in finalizer
    assert 'validate_resume(saved["fingerprint"]["digest"], fingerprint)' not in finalizer
    assert "if completion_path.exists() and completion_rebind is None:" in finalizer


def test_mapping_streams_preprocessing_and_applies_required_filters() -> None:
    inputs = MappingInputs(
        sample_id="sample-1",
        r1_paths=("reads_R1.fastq.gz",),
        r2_paths=("reads_R2.fastq.gz",),
        reference="combined.fa",
        output_bam="sample-1.bam",
        fastp_json="sample-1.fastp.json",
        fastp_html="sample-1.fastp.html",
    )

    command = build_mapping_command(inputs, threads=4)

    assert "fastp" in command
    assert "--qualified_quality_phred 20" in command
    assert "--unqualified_percent_limit 40" in command
    assert "--length_required 50" in command
    assert "--stdout" in command
    assert "bwa mem -p" in command
    assert "@RG\\tID:sample-1\\tSM:sample-1" in command
    assert "samtools view -F 3844 -q 20" in command
    assert "samtools sort -T sort_tmp/sample-1.name -n -@ 2 -u -" in command
    assert "samtools sort -T sort_tmp/sample-1.coordinate -@ 2 -u -" in command
    assert " -Ou " not in command
    assert "samtools fixmate -m" in command
    assert "samtools markdup -r" in command
    assert command.endswith("samtools index sample-1.bam")


def test_mapping_rejects_unbalanced_lanes() -> None:
    inputs = MappingInputs(
        sample_id="sample-1",
        r1_paths=("lane1_R1.fastq.gz", "lane2_R1.fastq.gz"),
        r2_paths=("lane1_R2.fastq.gz",),
        reference="combined.fa",
        output_bam="sample-1.bam",
        fastp_json="sample-1.fastp.json",
        fastp_html="sample-1.fastp.html",
    )

    try:
        build_mapping_command(inputs, threads=4)
    except ValueError as error:
        assert "balanced" in str(error)
    else:
        raise AssertionError("unbalanced lane inputs were accepted")
