from pathlib import Path

from dudleya_supplement.provenance import build_fingerprint, sha256_file


def _descendant_chain(inputs: dict[str, str]) -> tuple[str, ...]:
    upstream: dict[str, str] = {}
    digests = []
    for stage in ("eligibility", "variants", "statistics", "figures", "reports"):
        state = build_fingerprint(stage, inputs if stage == "eligibility" else {}, upstream, [stage], "commit")
        digests.append(state.digest)
        upstream = {stage: state.digest}
    return tuple(digests)


def test_miniature_dependency_chain_invalidates_every_descendant() -> None:
    fixture = Path(__file__).parent / "fixtures/mini"
    baseline = {
        "alignment": sha256_file(fixture / "alignment.fa"),
        "metadata": sha256_file(fixture / "samples.tsv"),
        "mask": sha256_file(fixture / "mask.bed"),
        "threshold": "dp5-gq20",
        "seed": "424200",
    }
    original = _descendant_chain(baseline)
    for key, changed in (
        ("alignment", "changed-canonical-input"),
        ("metadata", "changed-metadata"),
        ("mask", "changed-mask"),
        ("threshold", "dp10-gq30"),
        ("seed", "424201"),
    ):
        modified = {**baseline, key: changed}
        assert all(left != right for left, right in zip(original, _descendant_chain(modified), strict=True))


def test_repository_has_no_supplementary_symlink_into_forbidden_zones() -> None:
    root = Path(__file__).resolve().parents[3]
    forbidden = [root / "canonical_publication", root / "archive_noncanonical", root / "source_data"]
    for path in (root / "supplementary_analysis").rglob("*"):
        if path.is_symlink():
            target = path.resolve()
            assert not any(target == base or target.is_relative_to(base) for base in forbidden)


def test_supplementary_sources_have_no_workstation_absolute_path() -> None:
    root = Path(__file__).resolve().parents[3]
    forbidden = "/" + "/".join(("home", "neil"))
    for subtree in ("pipeline", "config", "README.md", "environment.yml", "run_pipeline.sh"):
        path = root / "supplementary_analysis" / subtree
        files = [path] if path.is_file() else list(path.rglob("*"))
        for file in files:
            if file.is_file() and "fixtures" not in file.parts:
                assert forbidden not in file.read_text(errors="ignore")
