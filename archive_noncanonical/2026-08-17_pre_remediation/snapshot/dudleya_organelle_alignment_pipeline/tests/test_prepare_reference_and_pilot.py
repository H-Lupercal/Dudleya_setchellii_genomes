import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.prepare_reference_and_pilot import (  # noqa: E402
    ReferenceValidationError,
    read_fasta_lengths,
    select_pilot_samples,
    validate_reference_records,
    write_pilot_samples,
)


def sample_row(
    sample_id,
    species,
    popcode,
    naming_profile="main_standard",
    analysis_status="include_primary_paired_end",
):
    return {
        "sample_id": sample_id,
        "batch": "QB3.Berkeley.251217/QB3.Results.260109",
        "naming_profile": naming_profile,
        "popcode": popcode,
        "species": species,
        "population_name": f"{popcode} population",
        "du_id": "Du-1",
        "lp_id": "LP_001",
        "sequencing_samples": "S1",
        "lanes": "L005",
        "r1_paths": f"{sample_id}_R1.fastq.gz",
        "r2_paths": f"{sample_id}_R2.fastq.gz",
        "r1_count": "1",
        "r2_count": "1",
        "pair_status": "complete",
        "metadata_status": "resolved",
        "analysis_status": analysis_status,
        "analysis_note": "",
    }


class ReferenceValidationTests(unittest.TestCase):
    def test_reads_fasta_lengths(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "reference.fa"
            fasta.write_text(">chloroplast\nAAAA\nAA\n>mitochondria\nCCC\n")

            lengths = read_fasta_lengths(fasta)

        self.assertEqual(lengths, {"chloroplast": 6, "mitochondria": 3})

    def test_validates_expected_reference_records(self):
        checks = validate_reference_records(
            {"chloroplast": 6, "mitochondria": 3},
            {"chloroplast": 6, "mitochondria": 3},
        )

        self.assertEqual([check.status for check in checks], ["PASS", "PASS"])

    def test_rejects_missing_reference_record(self):
        with self.assertRaises(ReferenceValidationError):
            validate_reference_records(
                {"chloroplast": 6},
                {"chloroplast": 6, "mitochondria": 3},
            )


class PilotSelectionTests(unittest.TestCase):
    def test_selects_diverse_complete_samples_and_skips_excluded_rows(self):
        rows = [
            sample_row("CY_RED_1", "D. cymosa", "CY_RED"),
            sample_row("CY_BAL_1", "D. cymosa", "CY_BAL"),
            sample_row("ABAB_MAD_1", "D. abramsii ssp. abramsii", "ABAB_MAD"),
            sample_row("BAI_1", "D. setchellii", "BAI"),
            sample_row(
                "DU-4A",
                "",
                "",
                naming_profile="initial_du_dash",
            ),
            sample_row(
                "DU014LP012",
                "",
                "",
                naming_profile="initial_du_lp",
            ),
            sample_row(
                "QUI1_missing",
                "D. setchellii",
                "QUI1",
                analysis_status="exclude_missing_mate",
            ),
        ]

        pilot = select_pilot_samples(rows, max_samples=6)

        selected = [row["sample_id"] for row in pilot]
        self.assertEqual(len(selected), 6)
        self.assertIn("DU-4A", selected)
        self.assertIn("DU014LP012", selected)
        self.assertIn("CY_BAL_1", selected)
        self.assertIn("ABAB_MAD_1", selected)
        self.assertIn("BAI_1", selected)
        self.assertNotIn("QUI1_missing", selected)
        self.assertTrue(all(row["pilot_reason"] for row in pilot))

    def test_write_pilot_samples_adds_reason_column(self):
        rows = [sample_row("CY_RED_1", "D. cymosa", "CY_RED")]
        pilot = select_pilot_samples(rows, max_samples=1)

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "pilot_samples.tsv"
            write_pilot_samples(out, pilot)
            with out.open() as handle:
                written = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(written[0]["sample_id"], "CY_RED_1")
        self.assertEqual(written[0]["pilot_reason"], "representative_D_cymosa")

    def test_fill_prefers_main_dataset_populations_after_initial_representatives(self):
        rows = [
            sample_row(
                "DU-4A",
                "",
                "",
                naming_profile="initial_du_dash",
            ),
            sample_row(
                "DU-173",
                "",
                "",
                naming_profile="initial_du_dash",
            ),
            sample_row(
                "DU014LP012",
                "",
                "",
                naming_profile="initial_du_lp",
            ),
            sample_row("CY_RED_1", "D. cymosa", "CY_RED"),
            sample_row("ABAB_MAD_1", "D. abramsii ssp. abramsii", "ABAB_MAD"),
            sample_row("BAI_1", "D. setchellii", "BAI"),
            sample_row("COM_1", "D. setchellii", "COM"),
        ]

        pilot = select_pilot_samples(rows, max_samples=6)

        selected = [row["sample_id"] for row in pilot]
        initial_du_dash_count = sum(
            1 for row in pilot if row["naming_profile"] == "initial_du_dash"
        )
        self.assertEqual(initial_du_dash_count, 1)
        self.assertIn("DU014LP012", selected)
        self.assertIn("COM_1", selected)

    def test_main_dataset_fill_round_robins_across_species_groups(self):
        rows = [
            sample_row(
                "DU-4A",
                "",
                "",
                naming_profile="initial_du_dash",
            ),
            sample_row(
                "DU014LP012",
                "",
                "",
                naming_profile="initial_du_lp",
            ),
            sample_row("CY_RED_1", "D. cymosa", "CY_RED"),
            sample_row("CY_BAL_1", "D. cymosa", "CY_BAL"),
            sample_row("CY_BGL_1", "D. cymosa", "CY_BGL"),
            sample_row("ABAB_MAD_1", "D. abramsii ssp. abramsii", "ABAB_MAD"),
            sample_row("ABBE_BERN_1", "D. abramsii ssp. bettinae", "ABBE_BERN"),
            sample_row("ABMU_HOR_1", "D. abramsii ssp. murina", "ABMU_HOR"),
            sample_row("BAI_1", "D. setchellii", "BAI"),
            sample_row("COM_1", "D. setchellii", "COM"),
        ]

        pilot = select_pilot_samples(rows, max_samples=8)

        selected = [row["sample_id"] for row in pilot]
        self.assertIn("CY_BAL_1", selected)
        self.assertIn("ABBE_BERN_1", selected)
        self.assertIn("COM_1", selected)


if __name__ == "__main__":
    unittest.main()
