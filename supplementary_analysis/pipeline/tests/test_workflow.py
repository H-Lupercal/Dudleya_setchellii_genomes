import inspect
import subprocess
import sys
from pathlib import Path

import pytest
from dudleya_supplement import stages, workflow
from dudleya_supplement.workflow import STAGES


def test_workflow_order_matches_approved_dependency_graph() -> None:
    assert STAGES == (
        "canonical_guard",
        "metadata",
        "identity",
        "sensitivity",
        "claims",
        "inheritance",
        "phase1_gate",
        "likelihood_mapping",
        "comparative_analyses",
        "figures",
        "reports",
        "acceptance",
        "canonical_guard_final",
    )


def test_dry_run_never_contains_preprocessing_or_mapping_commands() -> None:
    root = Path(__file__).resolve().parents[3]
    result = subprocess.run(
        [
            str(root / "supplementary_analysis/run_pipeline.sh"),
            "--config",
            "supplementary_analysis/config/supplementary_config.toml",
            "--run-id",
            "supplement-20260824",
            "--dry-run",
        ],
        cwd=root,
        capture_output=True,
        text=True,
        check=True,
    )
    output = result.stdout.lower()
    assert "stage canonical_guard" in output
    assert "stage canonical_guard_final" in output
    assert "fastp" not in output
    assert "bwa mem" not in output
    assert "map_samples.py" not in output


def test_every_declared_stage_has_an_implementation() -> None:
    source = inspect.getsource(stages.run_stage)
    assert "run_unimplemented" not in source


def test_from_stage_executes_upstream_resume_validation(monkeypatch) -> None:
    root = Path(__file__).resolve().parents[3]
    commands: list[list[str]] = []

    def record(command, **_kwargs):
        commands.append(command)

    monkeypatch.setattr(workflow.subprocess, "run", record)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            "--config",
            "supplementary_analysis/config/supplementary_config.toml",
            "--run-id",
            "supplement-20260824",
            "--resume",
            "--from-stage",
            "reports",
        ],
    )

    assert workflow.main() == 0
    assert len(commands) == len(STAGES)
    assert [command[command.index("--stage") + 1] for command in commands] == list(STAGES)
    assert all("--resume" in command for command in commands)
    assert root == Path.cwd()


def test_v26_config_rejects_a_different_run_id(monkeypatch) -> None:
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_pipeline.py",
            "--config",
            "supplementary_analysis/config/supplementary_config.v2.6.toml",
            "--run-id",
            "wrong-run",
            "--dry-run",
        ],
    )
    with pytest.raises(SystemExit, match="requires --run-id supplement-20260824-v26"):
        workflow.main()
