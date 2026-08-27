import csv
import hashlib
import importlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest
from organelle_pipeline import popgen


def test_pairwise_distance_counts_only_jointly_callable_substitutions() -> None:
    result = popgen.pairwise_sequence_distance("AaNC?T", "ATGCGN")

    assert result.differences == 1
    assert result.sites_compared == 3
    assert result.p_distance == pytest.approx(1 / 3)


def test_pairwise_distance_handles_identical_and_wholly_uncallable_sequences() -> None:
    identical = popgen.pairwise_sequence_distance("acgt", "ACGT")
    missing = popgen.pairwise_sequence_distance("NN?", "---")

    assert identical.differences == 0
    assert identical.sites_compared == 4
    assert identical.p_distance == 0.0
    assert missing.differences == 0
    assert missing.sites_compared == 0
    assert math.isnan(missing.p_distance)


def test_pairwise_distance_rejects_unequal_sequence_lengths() -> None:
    with pytest.raises(ValueError, match="equal lengths"):
        popgen.pairwise_sequence_distance("AC", "A")


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def _consensus_state(root: Path, run_id: str, organelle: str, outputs: tuple[Path, ...]) -> None:
    path = root / f"canonical_publication/provenance/runs/{run_id}/consensus/{organelle}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "status": "complete",
                "fingerprint": {"digest": f"consensus-{organelle}"},
                "outputs": {
                    output.relative_to(root).as_posix(): hashlib.sha256(output.read_bytes()).hexdigest() for output in outputs
                },
            }
        )
    )


def _read_tsv(path: Path) -> list[list[str]]:
    with path.open(newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def test_distance_stage_writes_ordered_symmetric_matrices_long_form_and_resume_state(tmp_path: Path) -> None:
    run_id = "miniature-run"
    metadata = (
        "sample_id\tpopcode\tspecies\n"
        "s2\tP1\tD. setchellii\n"
        "s1\tP1\tD. setchellii\n"
        "s3\tP2\tD. cymosa\n"
    )
    alignment = ">s1\nACGTN\n>s2\nATGTN\n>s3\nNNGTA\n"
    for organelle in ("chloroplast", "mitochondria"):
        metadata_path = _write(
            tmp_path / f"canonical_publication/metadata/qc/{run_id}/{organelle}_samples.tsv",
            metadata,
        )
        alignment_path = _write(
            tmp_path / f"canonical_publication/results/alignments/{run_id}/{organelle}.callable_alignment.fa",
            alignment,
        )
        _consensus_state(tmp_path, run_id, organelle, (alignment_path,))

    repository_root = Path(__file__).resolve().parents[3]
    script = repository_root / "canonical_publication/pipeline/scripts/sample_distances.py"
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(repository_root / "canonical_publication/pipeline/src")
    command = [
        sys.executable,
        str(script),
        "--repository-root",
        str(tmp_path),
        "--run-id",
        run_id,
    ]
    subprocess.run(command, cwd=repository_root, env=environment, check=True)

    for organelle in ("chloroplast", "mitochondria"):
        output = tmp_path / f"canonical_publication/results/supplement/{run_id}/pairwise_distances"
        differences = _read_tsv(output / f"{organelle}.sample_pairwise_differences.tsv")
        callable_sites = _read_tsv(output / f"{organelle}.sample_pairwise_callable_sites.tsv")
        long_form = _read_tsv(output / f"{organelle}.sample_pairwise_distances.tsv")

        assert differences == [
            ["sample_id", "s2", "s1", "s3"],
            ["s2", "0", "1", "0"],
            ["s1", "1", "0", "0"],
            ["s3", "0", "0", "0"],
        ]
        assert callable_sites == [
            ["sample_id", "s2", "s1", "s3"],
            ["s2", "4", "4", "2"],
            ["s1", "4", "4", "2"],
            ["s3", "2", "2", "3"],
        ]
        assert long_form[0] == ["organelle", "sample_1", "sample_2", "differences", "sites_compared", "p_distance"]
        assert long_form[1:] == [
            [organelle, "s2", "s1", "1", "4", "0.25"],
            [organelle, "s2", "s3", "0", "2", "0"],
            [organelle, "s1", "s3", "0", "2", "0"],
        ]
        state = json.loads(
            (tmp_path / f"canonical_publication/provenance/runs/{run_id}/distances/{organelle}.json").read_text()
        )
        assert state["status"] == "complete"
        assert state["sample_count"] == 3
        assert state["pairwise_comparison_count"] == 3
        assert len(state["outputs"]) == 3

    resumed = subprocess.run(
        [*command, "--resume"],
        cwd=repository_root,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
    assert "resume-valid sample distances chloroplast" in resumed.stdout
    assert "resume-valid sample distances mitochondria" in resumed.stdout


def test_distance_output_validator_checks_contract_and_summarizes_pairwise_counts(tmp_path: Path) -> None:
    distance_dir = tmp_path / "distances"
    differences = _write(
        distance_dir / "chloroplast.sample_pairwise_differences.tsv",
        "sample_id\ts1\ts2\ts3\ns1\t0\t1\t0\ns2\t1\t0\t2\ns3\t0\t2\t0\n",
    )
    callable_sites = _write(
        distance_dir / "chloroplast.sample_pairwise_callable_sites.tsv",
        "sample_id\ts1\ts2\ts3\ns1\t4\t4\t2\ns2\t4\t4\t2\ns3\t2\t2\t3\n",
    )
    long_form = _write(
        distance_dir / "chloroplast.sample_pairwise_distances.tsv",
        "organelle\tsample_1\tsample_2\tdifferences\tsites_compared\tp_distance\n"
        "chloroplast\ts1\ts2\t1\t4\t0.25\n"
        "chloroplast\ts1\ts3\t0\t2\t0\n"
        "chloroplast\ts2\ts3\t2\t2\t1\n",
    )
    distances = importlib.import_module("organelle_pipeline.distances")

    summary = distances.validate_pairwise_distance_outputs(
        "chloroplast",
        ("s1", "s2", "s3"),
        differences,
        callable_sites,
        long_form,
    )

    assert summary.sample_count == 3
    assert summary.pair_count == 3
    assert summary.minimum_differences == 0
    assert summary.median_differences == 1
    assert summary.maximum_differences == 2


def test_distance_output_validator_rejects_asymmetric_matrix(tmp_path: Path) -> None:
    distance_dir = tmp_path / "distances"
    differences = _write(
        distance_dir / "mitochondria.sample_pairwise_differences.tsv",
        "sample_id\ts1\ts2\ns1\t0\t1\ns2\t2\t0\n",
    )
    callable_sites = _write(
        distance_dir / "mitochondria.sample_pairwise_callable_sites.tsv",
        "sample_id\ts1\ts2\ns1\t4\t4\ns2\t4\t4\n",
    )
    long_form = _write(
        distance_dir / "mitochondria.sample_pairwise_distances.tsv",
        "organelle\tsample_1\tsample_2\tdifferences\tsites_compared\tp_distance\n"
        "mitochondria\ts1\ts2\t1\t4\t0.25\n",
    )
    distances = importlib.import_module("organelle_pipeline.distances")

    with pytest.raises(ValueError, match="symmetric"):
        distances.validate_pairwise_distance_outputs(
            "mitochondria",
            ("s1", "s2"),
            differences,
            callable_sites,
            long_form,
        )
