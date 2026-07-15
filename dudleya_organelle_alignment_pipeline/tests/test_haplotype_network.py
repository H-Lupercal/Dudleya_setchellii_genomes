import csv
import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.haplotype_network import (
    HaplotypeNetworkError,
    filter_complete_case_sites,
    validate_sample_metadata,
    write_network_input_fasta,
    write_network_metadata,
    write_network_site_table,
)


class CompleteCaseTests(unittest.TestCase):
    def test_filter_complete_case_sites_drops_any_non_acgt_column(self):
        records = [("S1", "ACNT"), ("S2", "ATGT"), ("S3", "ACGT")]

        filtered, kept, dropped = filter_complete_case_sites(records)

        self.assertEqual(
            filtered,
            [("S1", "ACT"), ("S2", "ATT"), ("S3", "ACT")],
        )
        self.assertEqual(kept, [0, 1, 3])
        self.assertEqual(dropped, [2])

    def test_validate_sample_metadata_rejects_mismatch(self):
        with self.assertRaisesRegex(HaplotypeNetworkError, "sample IDs"):
            validate_sample_metadata(["S1", "S2"], {"S1": {}, "S3": {}})

    def test_input_writers_record_site_filter_and_unresolved_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            fasta_path = output_dir / "network.fa"
            site_path = output_dir / "sites.tsv"
            metadata_path = output_dir / "metadata.tsv"
            records = [("S1", "ACT"), ("S2", "ATT")]
            source_rows = [
                {"organelle": "cpDNA", "position": str(position)}
                for position in (10, 20, 30, 40)
            ]
            metadata = {
                "S1": {
                    "species": "D. cymosa",
                    "popcode": "CY_RED",
                    "population_name": "Red",
                    "naming_profile": "main_standard",
                },
                "S2": {
                    "species": "",
                    "popcode": "",
                    "population_name": "",
                    "naming_profile": "initial_du_dash",
                },
            }

            write_network_input_fasta(fasta_path, records)
            write_network_site_table(site_path, source_rows, [0, 1, 3])
            write_network_metadata(metadata_path, records, metadata)

            fasta_text = fasta_path.read_text()
            with site_path.open(newline="") as handle:
                site_rows = list(csv.DictReader(handle, delimiter="\t"))
            with metadata_path.open(newline="") as handle:
                metadata_rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(
            [row["network_status"] for row in site_rows],
            ["retained", "retained", "dropped_missing", "retained"],
        )
        self.assertEqual(fasta_text, ">S1\nACT\n>S2\nATT\n")
        self.assertEqual(metadata_rows[1]["species_group"], "unresolved")


if __name__ == "__main__":
    unittest.main()
