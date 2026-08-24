from pathlib import Path

from dudleya_supplement.likelihood import build_likelihood_command, parse_likelihood_report
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
    assert parse_identical_sequence_map(text) == {
        "DU-229": "BAI_LP_110_Du-227",
        "CROB_LP_321_Du-396": "CROB_LP_173_Du-385",
    }


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
