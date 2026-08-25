import sys
from pathlib import Path

import pytest
from dudleya_supplement.likelihood import (
    build_likelihood_command,
    mask_restricted_sequences,
    parse_likelihood_diagnostics,
    parse_likelihood_report,
    run_command_logged,
)
from dudleya_supplement.phylogeny import (
    likelihood_decision,
    parse_identical_sequence_map,
    parse_split_nexus,
    supported_incompatible_pair,
)


def test_split_parser_and_incompatibility_use_bootstrap_frequencies(tmp_path: Path) -> None:
    nexus = tmp_path / "splits.nex"
    nexus.write_text(
        """#nexus
BEGIN Taxa;
DIMENSIONS ntax=4;
TAXLABELS
[1] 'A'
[2] 'B'
[3] 'C'
[4] 'D'
;
END;
BEGIN Splits;
DIMENSIONS ntax=4 nsplits=2;
FORMAT labels=no weights=yes;
MATRIX
30 1 2,
25 1 3,
;
END;
"""
    )
    taxa, splits = parse_split_nexus(nexus)
    assert taxa == ("A", "B", "C", "D")
    assert supported_incompatible_pair(taxa, splits, minimum_frequency=0.20) is not None


def test_neighbornet_requires_both_side_signal_and_supported_conflict() -> None:
    assert likelihood_decision(center_fraction=0.10, side_fraction=0.21, has_supported_conflict=True) == "RUN_NEIGHBORNET"
    assert likelihood_decision(center_fraction=0.16, side_fraction=0.30, has_supported_conflict=True) == "INSUFFICIENT_INFORMATION"
    assert likelihood_decision(center_fraction=0.10, side_fraction=0.21, has_supported_conflict=False) == "TREE_LIKE_NO_NETWORK"


def test_identical_tip_mapping_is_parsed_from_iqtree_notes() -> None:
    text = """NOTE: DU-229 is identical to BAI_LP_110_Du-227 but kept for subsequent analysis
NOTE: CROB_LP_321_Du-396 (identical to CROB_LP_173_Du-385) is ignored but added at the end
"""
    assert parse_identical_sequence_map(text) == {"CROB_LP_321_Du-396": "CROB_LP_173_Du-385"}


def test_likelihood_report_parses_all_seven_regions(tmp_path: Path) -> None:
    report = tmp_path / "test.iqtree"
    report.write_text(
        "Quartet support of areas 1-7 (mainly for clustered analysis):\n"
        " 100000  10000 ( 10.00) 20000 ( 20.00) 30000 ( 30.00) "
        "10000 ( 10.00) 10000 ( 10.00) 10000 ( 10.00) 10000 ( 10.00)\n"
        "Quartet resolution per sequence\n"
        "Overall quartet resolution:\n"
        "Number of fully resolved  quartets (regions 1+2+3): 60000 (=60.00%)\n"
        "Number of partly resolved quartets (regions 4+5+6): 30000 (=30.00%)\n"
        "Number of unresolved      quartets (region 7)     : 10000 (=10.00%)\n"
    )
    values = parse_likelihood_report(report)
    assert values["region_7_fraction"] == 0.10
    assert values["side_fraction"] == 0.30


def test_likelihood_report_accepts_iqtree3_percentage_spacing(tmp_path: Path) -> None:
    report = tmp_path / "iqtree3.iqtree"
    report.write_text(
        "Quartet support of areas 1-7 (mainly for clustered analysis):\n"
        "    100000   32786   (32.79 ) 32658   (32.66 ) 32670   (32.67 ) "
        "202 (0.20  ) 190 (0.19  ) 196 (0.20  ) 1298 (1.30  )\n"
        "Quartet resolution per sequence\n"
        "Overall quartet resolution:\n"
        "Number of fully resolved  quartets (regions 1+2+3): 98114 (=98.11%)\n"
        "Number of partly resolved quartets (regions 4+5+6): 588 (=0.59%)\n"
        "Number of unresolved      quartets (region 7)     : 1298 (=1.30%)\n"
    )
    values = parse_likelihood_report(report)
    assert values["region_1_count"] == 32786
    assert values["region_7_fraction"] == pytest.approx(0.013)


def test_likelihood_mapping_declares_dna_sequence_type() -> None:
    command = build_likelihood_command(
        alignment=Path("canonical/alignment.fa"),
        model="TVM+F+I+R4",
        quartets=100_000,
        seed=271828,
        prefix=Path("supplement/work/chloroplast"),
        threads=8,
    )
    assert command[command.index("-st") + 1] == "DNA"


def test_likelihood_command_output_is_redirected_to_work_log(tmp_path: Path) -> None:
    log = tmp_path / "screen.log"
    run_command_logged([sys.executable, "-c", "print('quartet diagnostics')"], cwd=tmp_path, log=log)
    assert log.read_text() == "quartet diagnostics\n"


def test_likelihood_diagnostics_capture_composition_and_ambiguity(tmp_path: Path) -> None:
    log = tmp_path / "mitochondria.log"
    log.write_text(
        "Alignment has 271 sequences with 243359 columns, 2405 distinct patterns\n"
        "WARNING: 271 sequences contain more than 50% gaps/ambiguity\n"
        "****  TOTAL 82.33%  269 sequences failed composition chi2 test (p-value<5%; df=3)\n"
    )
    assert parse_likelihood_diagnostics(log) == {
        "alignment_sequence_count": 271,
        "composition_failed_count": 269,
        "over_50pct_ambiguity_count": 271,
    }


def test_mask_restricted_alignment_uses_exact_zero_based_half_open_intervals() -> None:
    restricted = mask_restricted_sequences(
        {"sample-a": "AACCGG", "sample-b": "TTGGCC"},
        [(1, 3), (4, 6)],
        expected_length=4,
    )
    assert restricted == {"sample-a": "ACGG", "sample-b": "TGCC"}
    assert len(restricted) == 2


def test_mitochondrial_sensitivity_uses_matched_primary_likelihood_settings() -> None:
    command = build_likelihood_command(
        alignment=Path("supplement/mitochondria.mask_restricted.fa"),
        model="TPM3u+F+I",
        quartets=100_000,
        seed=314159,
        prefix=Path("supplement/work/mitochondria.mask_restricted"),
        threads=8,
    )
    assert command[command.index("-m") + 1] == "TPM3u+F+I"
    assert command[command.index("-lmap") + 1] == "100000"
    assert command[command.index("-seed") + 1] == "314159"


def test_mask_restricted_alignment_rejects_wrong_length() -> None:
    with pytest.raises(ValueError, match="expected 5"):
        mask_restricted_sequences({"sample": "AACCGG"}, [(1, 3), (4, 6)], expected_length=5)
