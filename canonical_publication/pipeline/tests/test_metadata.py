from pathlib import Path

from organelle_pipeline.metadata import discover_samples, read_population_codes


def test_manifest_is_regenerated_from_fastq_pairs_and_population_codes(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    for mate in (1, 2):
        (raw / f"POP_LP_001_Du-1_S1_L001_R{mate}_001.fastq.gz").touch()
    (raw / "ORPHAN_LP_002_Du-2_S2_L001_R1_001.fastq.gz").touch()
    codes = tmp_path / "codes.csv"
    codes.write_text("Species,Population Name,Code\nD. test,Test Population,POP\n")

    populations = read_population_codes(codes)
    samples = discover_samples(raw, populations)

    assert len(samples) == 2
    by_id = {sample.sample_id: sample for sample in samples}
    assert by_id["POP_LP_001_Du-1"].pair_status == "complete"
    assert by_id["POP_LP_001_Du-1"].popcode == "POP"
    assert by_id["ORPHAN_LP_002_Du-2"].pair_status == "missing_R2"


def test_multiple_lanes_are_retained_as_ordered_paths(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    for lane in ("L001", "L002"):
        for mate in (1, 2):
            (raw / f"POP_LP_001_Du-1_S1_{lane}_R{mate}_001.fastq.gz").touch()

    sample = discover_samples(raw, {})[0]

    assert len(sample.r1_paths) == 2
    assert len(sample.r2_paths) == 2
    assert sample.pair_status == "complete"


def test_du_numbered_ids_use_source_tables_declared_duse_default(tmp_path: Path) -> None:
    from organelle_pipeline.metadata import PopulationCode

    populations = {"DUSE": PopulationCode("DUSE", "D. setchellii", "default population")}
    raw = tmp_path / "raw"
    raw.mkdir()
    for sample in ("DU-173", "DU014LP012"):
        for mate in (1, 2):
            (raw / f"{sample}_S1_L001_R{mate}_001.fastq.gz").touch()
    observed = {sample.sample_id: sample.popcode for sample in discover_samples(raw, populations)}
    assert observed == {"DU-173": "DUSE", "DU014LP012": "DUSE"}


def test_population_table_header_materializes_declared_duse_default(tmp_path: Path) -> None:
    codes = tmp_path / "codes.csv"
    codes.write_text("Species,Population Name,Code (if it doesn't start with a TWO letter code = DUSE)\n,Boulder Ridge,CY_BOU\n")
    populations = read_population_codes(codes)
    assert populations["DUSE"].species == "D. setchellii"
    assert populations["CY_BOU"].species == "D. cymosa"


def test_duplicate_population_code_preserves_all_conflicting_source_labels(tmp_path: Path) -> None:
    codes = tmp_path / "codes.csv"
    codes.write_text("Species,Population Name,Code\nD. setchellii,First label,TUL2\nD. setchellii,Second label,TUL2\n")
    assert read_population_codes(codes)["TUL2"].population_name == "First label | Second label"
