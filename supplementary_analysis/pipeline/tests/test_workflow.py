import inspect
import subprocess
from pathlib import Path

from dudleya_supplement import stages
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
