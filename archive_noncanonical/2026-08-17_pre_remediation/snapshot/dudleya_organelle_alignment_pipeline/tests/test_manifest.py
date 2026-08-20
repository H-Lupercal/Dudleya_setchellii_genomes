import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.manifest import (  # noqa: E402
    build_manifest,
    discover_fastq_paths,
    load_population_codes,
    parse_fastq_path,
    write_manifest_outputs,
)


class FastqNameParsingTests(unittest.TestCase):
    def test_parses_initial_du_dash_names(self):
        record = parse_fastq_path(
            Path(
                "genomicsDrive_data_dump/QB3.Berkeley.241122/"
                "QB3.Dudleya.Results.250118/DU-4A_S68_L008_R1_001.fastq.gz"
            )
        )

        self.assertEqual(record.sample_id, "DU-4A")
        self.assertEqual(record.naming_profile, "initial_du_dash")
        self.assertEqual(record.read, "R1")
        self.assertEqual(record.sequencing_sample, "S68")
        self.assertEqual(record.lane, "L008")
        self.assertEqual(record.du_id, "DU-4A")
        self.assertEqual(record.popcode, "")

    def test_parses_initial_du_lp_names(self):
        record = parse_fastq_path(
            Path(
                "genomicsDrive_data_dump/QB3.Berkeley.250811/"
                "QB3.250916.Results5genomes/DU014LP012_S4_L005_R2_001.fastq.gz"
            )
        )

        self.assertEqual(record.sample_id, "DU014LP012")
        self.assertEqual(record.naming_profile, "initial_du_lp")
        self.assertEqual(record.read, "R2")
        self.assertEqual(record.du_id, "DU014")
        self.assertEqual(record.lp_id, "LP012")
        self.assertEqual(record.popcode, "")

    def test_parses_main_standard_names(self):
        record = parse_fastq_path(
            Path(
                "genomicsDrive_data_dump/QB3.Berkeley.251217/QB3.Results.260109/"
                "CY_RED/CY_RED_LP_202_Du-561_S192_L005_R1_001.fastq.gz"
            )
        )

        self.assertEqual(record.sample_id, "CY_RED_LP_202_Du-561")
        self.assertEqual(record.naming_profile, "main_standard")
        self.assertEqual(record.read, "R1")
        self.assertEqual(record.popcode, "CY_RED")
        self.assertEqual(record.lp_id, "LP_202")
        self.assertEqual(record.du_id, "Du-561")

    def test_parses_main_standard_names_with_lp_du_hyphen_separator(self):
        record = parse_fastq_path(
            Path(
                "genomicsDrive_data_dump/QB3.Berkeley.251217/QB3.Results.260109/"
                "CY_ALA/CY_ALA_LP_298-Du-767_S289_L005_R1_001.fastq.gz"
            )
        )

        self.assertEqual(record.sample_id, "CY_ALA_LP_298-Du-767")
        self.assertEqual(record.naming_profile, "main_standard")
        self.assertEqual(record.popcode, "CY_ALA")
        self.assertEqual(record.lp_id, "LP_298")
        self.assertEqual(record.du_id, "Du-767")


class ManifestTests(unittest.TestCase):
    def test_discovery_finds_chunked_fastq_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            fastq_path = (
                Path(tmp)
                / "QB3.Berkeley.250811"
                / "QB3.250916.Results5genomes"
                / "DU014LP012_S4_L005_R1_001.fastq-011.gz"
            )
            fastq_path.parent.mkdir(parents=True)
            fastq_path.write_text("")

            discovered = discover_fastq_paths(Path(tmp))

        self.assertEqual(discovered, [fastq_path])

    def test_population_codes_are_loaded_by_code(self):
        with tempfile.TemporaryDirectory() as tmp:
            csv_path = Path(tmp) / "population_codes.csv"
            with csv_path.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "Species",
                        "Population Name",
                        "Code (if it doesn't start with a TWO letter code = DUSE)",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "Species": "D. cymosa",
                        "Population Name": "Red Mountain",
                        "Code (if it doesn't start with a TWO letter code = DUSE)": "CY_RED",
                    }
                )

            codes = load_population_codes(csv_path)

        self.assertIn("CY_RED", codes)
        self.assertEqual(codes["CY_RED"].species, "D. cymosa")
        self.assertEqual(codes["CY_RED"].population_name, "Red Mountain")

    def test_build_manifest_pairs_reads_and_marks_unresolved_metadata(self):
        paths = [
            Path(
                "genomicsDrive_data_dump/QB3.Berkeley.251217/QB3.Results.260109/"
                "CY_RED/CY_RED_LP_202_Du-561_S192_L005_R1_001.fastq.gz"
            ),
            Path(
                "genomicsDrive_data_dump/QB3.Berkeley.251217/QB3.Results.260109/"
                "CY_RED/CY_RED_LP_202_Du-561_S192_L005_R2_001.fastq.gz"
            ),
            Path(
                "genomicsDrive_data_dump/QB3.Berkeley.241122/"
                "QB3.Dudleya.Results.250118/DU-4A_S68_L008_R1_001.fastq.gz"
            ),
        ]
        population_codes = {
            "CY_RED": type(
                "PopulationCode",
                (),
                {
                    "code": "CY_RED",
                    "species": "D. cymosa",
                    "population_name": "Red Mountain",
                },
            )()
        }

        rows, issues = build_manifest(paths, population_codes)

        by_sample = {row.sample_id: row for row in rows}
        self.assertEqual(by_sample["CY_RED_LP_202_Du-561"].pair_status, "complete")
        self.assertEqual(by_sample["CY_RED_LP_202_Du-561"].species, "D. cymosa")
        self.assertEqual(by_sample["CY_RED_LP_202_Du-561"].metadata_status, "resolved")
        self.assertEqual(
            by_sample["CY_RED_LP_202_Du-561"].analysis_status,
            "include_primary_paired_end",
        )
        self.assertEqual(by_sample["DU-4A"].pair_status, "missing_R2")
        self.assertEqual(by_sample["DU-4A"].metadata_status, "unresolved_initial_sample")
        self.assertEqual(by_sample["DU-4A"].analysis_status, "exclude_missing_mate")
        self.assertIn("single-end", by_sample["DU-4A"].analysis_note)
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].issue_type, "missing_R2")

    def test_write_manifest_outputs_splits_primary_and_excluded_samples(self):
        paths = [
            Path(
                "genomicsDrive_data_dump/QB3.Berkeley.251217/QB3.Results.260109/"
                "CY_RED/CY_RED_LP_202_Du-561_S192_L005_R1_001.fastq.gz"
            ),
            Path(
                "genomicsDrive_data_dump/QB3.Berkeley.251217/QB3.Results.260109/"
                "CY_RED/CY_RED_LP_202_Du-561_S192_L005_R2_001.fastq.gz"
            ),
            Path(
                "genomicsDrive_data_dump/QB3.Berkeley.251217/QB3.Results.260109/"
                "QUI1/QUI1_LP_256_Du-655_S264_L005_R1_001.fastq.gz"
            ),
        ]

        rows, issues = build_manifest(paths, {})

        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            write_manifest_outputs(rows, issues, out)
            with (out / "analysis_samples.tsv").open() as handle:
                analysis_rows = list(csv.DictReader(handle, delimiter="\t"))
            with (out / "excluded_samples.tsv").open() as handle:
                excluded_rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual([row["sample_id"] for row in analysis_rows], ["CY_RED_LP_202_Du-561"])
        self.assertEqual([row["sample_id"] for row in excluded_rows], ["QUI1_LP_256_Du-655"])
        self.assertEqual(excluded_rows[0]["analysis_status"], "exclude_missing_mate")
        self.assertIn("single-end", excluded_rows[0]["analysis_note"])


if __name__ == "__main__":
    unittest.main()
