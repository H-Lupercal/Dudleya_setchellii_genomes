from pathlib import Path

import pytest
from organelle_pipeline.logs import portable_command_log
from organelle_pipeline.provenance import (
    StaleOutputError,
    build_stage_fingerprint,
    build_stage_fingerprint_from_hashes,
    mapping_pipeline_code_digest,
    pipeline_code_digest,
    runtime_pipeline_code_digest,
    runtime_provenance,
    sha256_json,
    validate_recorded_tool_versions,
    validate_resume,
    validate_saved_outputs,
)


def test_input_content_change_invalidates_stage_and_descendant(tmp_path: Path) -> None:
    source = tmp_path / "source.txt"
    source.write_text("first")
    upstream = build_stage_fingerprint("reference", [source], {}, [])
    child = build_stage_fingerprint("mapping", [source], {"reference": upstream.digest}, [])

    source.write_text("second")
    changed = build_stage_fingerprint("reference", [source], {}, [])
    changed_child = build_stage_fingerprint("mapping", [source], {"reference": changed.digest}, [])

    assert changed.digest != upstream.digest
    assert changed_child.digest != child.digest


def test_resume_rejects_stale_fingerprint(tmp_path: Path) -> None:
    source = tmp_path / "source.txt"
    source.write_text("current")
    current = build_stage_fingerprint("qc", [source], {}, ["samtools depth"])

    with pytest.raises(StaleOutputError, match="stale"):
        validate_resume("not-the-current-digest", current)


def test_resume_accepts_exact_fingerprint(tmp_path: Path) -> None:
    source = tmp_path / "source.txt"
    source.write_text("current")
    current = build_stage_fingerprint("qc", [source], {}, ["samtools depth"])
    validate_resume(current.digest, current)


def test_precomputed_source_hash_change_invalidates_without_rereading_files() -> None:
    first = build_stage_fingerprint_from_hashes("mapping", {"reads.fastq.gz": "aaa"}, {}, ["map"])
    changed = build_stage_fingerprint_from_hashes("mapping", {"reads.fastq.gz": "bbb"}, {}, ["map"])
    assert first.digest != changed.digest


def test_structured_digest_is_order_independent_but_content_sensitive() -> None:
    assert sha256_json({"a": 1, "b": 2}) == sha256_json({"b": 2, "a": 1})
    assert sha256_json({"a": 1}) != sha256_json({"a": 2})


def test_portable_command_log_replaces_repository_root_and_rejects_other_workstation_paths(tmp_path: Path) -> None:
    repository_root = tmp_path / "repo"
    repository_root.mkdir()

    normalized = portable_command_log(
        f"Working directory: {repository_root}\n--out {repository_root}/canonical_publication/work/result\n",
        repository_root,
    )

    assert str(repository_root) not in normalized
    assert normalized == ("Working directory: ${REPOSITORY_ROOT}\n--out ${REPOSITORY_ROOT}/canonical_publication/work/result\n")
    uncontrolled_path = "/" + "tmp/uncontrolled-tool-cache"
    with pytest.raises(ValueError, match="absolute workstation path"):
        portable_command_log(f"cache: {uncontrolled_path}\n", repository_root)


def test_pipeline_code_digest_changes_when_implementation_changes(tmp_path: Path) -> None:
    source = tmp_path / "canonical_publication" / "pipeline" / "src"
    source.mkdir(parents=True)
    module = source / "module.py"
    module.write_text("VALUE = 1\n")
    first = pipeline_code_digest(tmp_path)
    module.write_text("VALUE = 2\n")
    assert pipeline_code_digest(tmp_path) != first


def test_runtime_code_digest_tracks_loaded_stage_code_not_unrelated_downstream_code(tmp_path: Path) -> None:
    pipeline = tmp_path / "canonical_publication" / "pipeline"
    pipeline.mkdir(parents=True)
    stage = pipeline / "stage.py"
    downstream = pipeline / "downstream.py"
    stage.write_text("VALUE = 1\n")
    downstream.write_text("VALUE = 1\n")

    first = runtime_pipeline_code_digest(tmp_path, loaded_files=(stage,))
    downstream.write_text("VALUE = 2\n")
    assert runtime_pipeline_code_digest(tmp_path, loaded_files=(stage,)) == first
    stage.write_text("VALUE = 2\n")
    assert runtime_pipeline_code_digest(tmp_path, loaded_files=(stage,)) != first


def test_mapping_code_digest_is_shared_by_executor_and_validator(tmp_path: Path) -> None:
    pipeline = tmp_path / "canonical_publication/pipeline"
    scripts = pipeline / "scripts"
    package = pipeline / "src/organelle_pipeline"
    scripts.mkdir(parents=True)
    package.mkdir(parents=True)
    relevant = (
        scripts / "map_samples.py",
        package / "__init__.py",
        package / "commands.py",
        package / "inventory.py",
        package / "metadata.py",
        package / "paths.py",
        package / "provenance.py",
    )
    for path in relevant:
        path.write_text(f"# {path.name}\n")
    unrelated = package / "consensus.py"
    unrelated.write_text("# old downstream code\n")

    first = mapping_pipeline_code_digest(tmp_path)
    unrelated.write_text("# new downstream code\n")
    assert mapping_pipeline_code_digest(tmp_path) == first
    (package / "commands.py").write_text("# changed mapping command\n")
    assert mapping_pipeline_code_digest(tmp_path) != first


def test_runtime_provenance_records_code_commit_and_executable_version(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    pipeline = tmp_path / "canonical_publication" / "pipeline"
    pipeline.mkdir(parents=True)
    (pipeline / "runner.py").write_text("pass\n")

    class Completed:
        returncode = 0
        stdout = "tool 1.2.3\n"
        stderr = ""

    monkeypatch.setattr("organelle_pipeline.provenance.subprocess.run", lambda *args, **kwargs: Completed())
    recorded = runtime_provenance(tmp_path, {"tool": ("tool", "--version")})
    assert recorded["provenance:git_commit"] == "tool 1.2.3"
    assert recorded["provenance:tool:tool"] == "tool 1.2.3"
    assert "provenance:executable_sha256:tool" in recorded
    assert recorded["provenance:pipeline_code"] == runtime_pipeline_code_digest(tmp_path)


def test_resume_rejects_modified_saved_output(tmp_path: Path) -> None:
    output = tmp_path / "result.tsv"
    output.write_text("original\n")
    saved = {"outputs": {"result.tsv": "not-the-current-checksum"}}
    with pytest.raises(StaleOutputError, match="checksum"):
        validate_saved_outputs(tmp_path, saved)


def test_tool_version_validation_rejects_changed_or_missing_required_tool() -> None:
    current = {"fastp": "fastp 1.3.6", "bwa": "Version: 0.7.19", "samtools": "samtools 1.23.1"}
    with pytest.raises(StaleOutputError, match="bwa"):
        validate_recorded_tool_versions(
            {"fastp": "fastp 1.3.6", "bwa": "Version: 0.7.18", "samtools": "samtools 1.23.1"},
            current,
        )
    with pytest.raises(StaleOutputError, match="samtools"):
        validate_recorded_tool_versions({"fastp": "fastp 1.3.6", "bwa": "Version: 0.7.19"}, current)


def test_tool_version_validation_can_transparently_reconcile_legacy_missing_bwa() -> None:
    current = {"fastp": "fastp 1.3.6", "bwa": "Version: 0.7.19", "samtools": "samtools 1.23.1"}
    resolved, notes = validate_recorded_tool_versions(
        {"fastp": "fastp 1.3.6", "bwa": "", "samtools": "samtools 1.23.1"},
        current,
        reconcilable_missing={"bwa"},
    )
    assert resolved == current
    assert notes == ("bwa version absent in provisional state; resolved from pinned executable during reconciliation",)


def test_miniature_dependency_chain_invalidates_for_reference_mask_threshold_and_sample_table_changes() -> None:
    def chain(reference: str, mask: str, threshold: str, samples: str) -> dict[str, str]:
        references = build_stage_fingerprint_from_hashes("references", {"reference.fa": reference}, {}, ["select"])
        mapping = build_stage_fingerprint_from_hashes(
            "mapping",
            {"reference.fa": reference, "reads.fastq": "reads"},
            {"references": references.digest},
            ["map"],
        )
        qc = build_stage_fingerprint_from_hashes(
            "qc",
            {"threshold": threshold, "samples.tsv": samples},
            {"mapping": mapping.digest},
            ["eligibility"],
        )
        variants = build_stage_fingerprint_from_hashes(
            "variants",
            {"mask.bed": mask},
            {"qc": qc.digest},
            ["call"],
        )
        report = build_stage_fingerprint_from_hashes(
            "report",
            {},
            {"variants": variants.digest},
            ["publish"],
        )
        return {
            stage: fingerprint.digest
            for stage, fingerprint in {
                "references": references,
                "mapping": mapping,
                "qc": qc,
                "variants": variants,
                "report": report,
            }.items()
        }

    baseline = chain("ref-v1", "mask-v1", "dp5=0.8", "samples-v1")
    assert chain("ref-v2", "mask-v1", "dp5=0.8", "samples-v1")["mapping"] != baseline["mapping"]
    assert chain("ref-v1", "mask-v2", "dp5=0.8", "samples-v1")["variants"] != baseline["variants"]
    assert chain("ref-v1", "mask-v1", "dp5=0.9", "samples-v1")["qc"] != baseline["qc"]
    assert chain("ref-v1", "mask-v1", "dp5=0.8", "samples-v2")["qc"] != baseline["qc"]
    for changed in (
        chain("ref-v2", "mask-v1", "dp5=0.8", "samples-v1"),
        chain("ref-v1", "mask-v2", "dp5=0.8", "samples-v1"),
        chain("ref-v1", "mask-v1", "dp5=0.9", "samples-v1"),
        chain("ref-v1", "mask-v1", "dp5=0.8", "samples-v2"),
    ):
        assert changed["report"] != baseline["report"]
