import csv
import subprocess
from pathlib import Path

import pytest
from organelle_pipeline.paths import (
    CanonicalPathError,
    assert_canonical_path,
    repository_relative,
    validate_run_id,
)


def test_archive_path_is_rejected(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    archive = repo / "archive_noncanonical" / "old.tsv"
    archive.parent.mkdir(parents=True)
    archive.touch()

    with pytest.raises(CanonicalPathError, match="noncanonical archive"):
        assert_canonical_path(archive, repo)


def test_symlink_into_archive_is_rejected(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    archive = repo / "archive_noncanonical" / "old.tsv"
    canonical = repo / "canonical_publication"
    archive.parent.mkdir(parents=True)
    canonical.mkdir(parents=True)
    archive.touch()
    link = canonical / "bad-link.tsv"
    link.symlink_to(archive)

    with pytest.raises(CanonicalPathError, match="noncanonical archive"):
        assert_canonical_path(link, repo)


def test_source_and_canonical_paths_are_allowed(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    source = repo / "source_data" / "reads.fastq.gz"
    result = repo / "canonical_publication" / "results" / "qc.tsv"

    assert assert_canonical_path(source, repo) == source.resolve()
    assert assert_canonical_path(result, repo) == result.resolve()


def test_repository_relative_accepts_relative_cli_paths(tmp_path: Path) -> None:
    repo = tmp_path / "repo"
    repo.mkdir()
    assert repository_relative("canonical_publication/config.toml", repo) == Path("canonical_publication/config.toml")


def test_run_id_cannot_escape_canonical_output_directories() -> None:
    assert validate_run_id("publication-20260817") == "publication-20260817"
    with pytest.raises(CanonicalPathError, match="run ID"):
        validate_run_id("../../archive_noncanonical")


def test_archive_manifest_represents_every_legacy_artifact() -> None:
    root = Path(__file__).resolve().parents[3]
    snapshot = root / "archive_noncanonical/2026-08-17_pre_remediation/snapshot"
    with (root / "archive_noncanonical/2026-08-17_pre_remediation/manifest.tsv").open(newline="") as handle:
        manifested = {row["archived_path"] for row in csv.DictReader(handle, delimiter="\t")}
    observed = {path.relative_to(root).as_posix() for path in snapshot.rglob("*") if path.is_file() or path.is_symlink()}
    assert observed == manifested


def test_no_canonical_symlink_resolves_into_archive() -> None:
    root = Path(__file__).resolve().parents[3]
    archive = (root / "archive_noncanonical").resolve()
    links = [path for path in (root / "canonical_publication").rglob("*") if path.is_symlink()]
    assert all(not path.resolve(strict=False).is_relative_to(archive) for path in links)


def test_canonical_text_has_no_absolute_workstation_path() -> None:
    root = Path(__file__).resolve().parents[3]
    suffixes = {".py", ".sh", ".toml", ".md", ".tsv", ".json", ".yml", ".yaml", ".txt", ".gff3", ".nex"}
    offenders = []
    for path in (root / "canonical_publication").rglob("*"):
        if not path.is_file() or path.suffix not in suffixes or "work" in path.parts:
            continue
        text = path.read_text(errors="ignore")
        absolute_prefixes = tuple("/" + value for value in ("home/", "Users/", "tmp/"))
        if any(prefix in text for prefix in absolute_prefixes):
            offenders.append(path.relative_to(root).as_posix())
    assert offenders == []


def test_only_audit_reporting_and_path_guard_name_the_archive() -> None:
    root = Path(__file__).resolve().parents[3]
    allowed = {"build_reports.py", "paths.py"}
    offenders = []
    for path in (root / "canonical_publication/pipeline").rglob("*.py"):
        if "archive_noncanonical/" in path.read_text() and path.name not in allowed and "tests" not in path.parts:
            offenders.append(path.relative_to(root).as_posix())
    assert offenders == []


def test_large_archive_matrices_are_local_only_but_small_reports_remain_visible() -> None:
    root = Path(__file__).resolve().parents[3]
    snapshot = root / "archive_noncanonical/2026-08-17_pre_remediation/snapshot/full_pipeline_run/results"
    large_derived = (
        snapshot / "06_all_sample_alignment/qc/DU079LP010.depth.tsv",
        snapshot / "11_callable_consensus/cpDNA.primary.callable_sites.tsv",
    )
    small_report = snapshot / "06_all_sample_alignment/all_sample_alignment_report.md"

    for path in large_derived:
        assert path.is_file()
        assert (
            subprocess.run(
                ["git", "check-ignore", "--quiet", str(path.relative_to(root))],
                cwd=root,
                check=False,
            ).returncode
            == 0
        )
    assert small_report.is_file()
    assert (
        subprocess.run(
            ["git", "check-ignore", "--quiet", str(small_report.relative_to(root))],
            cwd=root,
            check=False,
        ).returncode
        == 1
    )


def test_canonical_reference_masks_are_not_hidden_by_plink_bed_ignore_rule() -> None:
    root = Path(__file__).resolve().parents[3]
    mask = root / "canonical_publication/references/masks/chloroplast_population_sites.bed"

    assert mask.is_file()
    assert (
        subprocess.run(
            ["git", "check-ignore", "--quiet", str(mask.relative_to(root))],
            cwd=root,
            check=False,
        ).returncode
        == 1
    )
