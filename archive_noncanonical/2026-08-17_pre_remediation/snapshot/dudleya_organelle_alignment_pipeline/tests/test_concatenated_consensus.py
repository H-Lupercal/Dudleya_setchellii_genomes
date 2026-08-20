import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.concatenated_consensus import (
    ConcatenatedConsensusError,
    concatenate_consensus_alignments,
    read_fasta_alignment,
    run_concatenation,
)


class ConcatenatedConsensusBuildTests(unittest.TestCase):
    def test_concatenates_by_sample_id_in_cpdna_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cpdna_path = root / "cpDNA.fa"
            cpdna_path.write_text(">S1\nACGN\n>S2\nTTAA\n")
            mtdna_path = root / "mtDNA.fa"
            mtdna_path.write_text(">S2\nGG\n>S1\nNC\n")

            alignment = concatenate_consensus_alignments(cpdna_path, mtdna_path)

        self.assertEqual(alignment.sample_names, ("S1", "S2"))
        self.assertEqual(
            alignment.sequences,
            {"S1": "ACGNNC", "S2": "TTAAGG"},
        )
        self.assertEqual(alignment.cpdna_length, 4)
        self.assertEqual(alignment.mtdna_length, 2)
        self.assertEqual(alignment.combined_length, 6)
        self.assertEqual(alignment.mtdna_start, 5)

    def test_rejects_duplicate_fasta_identifiers(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cpdna_path = root / "cpDNA.fa"
            cpdna_path.write_text(">S1\nAAAA\n>S1\nCCCC\n")
            mtdna_path = root / "mtDNA.fa"
            mtdna_path.write_text(">S1\nGG\n")

            with self.assertRaisesRegex(
                ConcatenatedConsensusError,
                "Duplicate FASTA identifier S1",
            ):
                concatenate_consensus_alignments(cpdna_path, mtdna_path)


class ConcatenatedConsensusOutputTests(unittest.TestCase):
    def test_run_writes_combined_fasta_compatibility_summary_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cpdna_path = root / "cpDNA.fa"
            cpdna_path.write_text(">S1\nACGN\n>S2\nTTAA\n")
            mtdna_path = root / "mtDNA.fa"
            mtdna_path.write_text(">S2\nGG\n>S1\nNC\n")
            output_dir = root / "output"

            result = run_concatenation(
                cpdna_path=cpdna_path,
                mtdna_path=mtdna_path,
                output_dir=output_dir,
                run_label="primary",
            )
            written = read_fasta_alignment(result.fasta_path)
            compatibility = (
                output_dir / "primary.callable_consensus_summary.tsv"
            ).read_text()
            detailed = (
                output_dir / "primary.concatenated_consensus_summary.tsv"
            ).read_text()
            report = (
                output_dir / "primary.concatenated_consensus_report.md"
            ).read_text()

        self.assertEqual(result.sample_count, 2)
        self.assertEqual(result.combined_length, 6)
        self.assertEqual(written.sample_names, ("S1", "S2"))
        self.assertEqual(written.sequences["S1"], "ACGNNC")
        self.assertIn("cpDNA_mtDNA\tcpdna_then_mtdna\t2\t6", compatibility)
        self.assertIn("cpDNA_end\tmtDNA_start", detailed)
        self.assertIn("\t4\t5\t", detailed)
        self.assertIn("cpDNA positions: 1-4", report)
        self.assertIn("mtDNA positions: 5-6", report)

    def test_rejects_mismatched_sample_identifier_sets(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cpdna_path = root / "cpDNA.fa"
            cpdna_path.write_text(">S1\nAAAA\n>S2\nCCCC\n")
            mtdna_path = root / "mtDNA.fa"
            mtdna_path.write_text(">S1\nGG\n>S3\nTT\n")

            with self.assertRaisesRegex(
                ConcatenatedConsensusError,
                "Sample identifier mismatch",
            ):
                concatenate_consensus_alignments(cpdna_path, mtdna_path)

    def test_rejects_inconsistent_record_lengths_in_either_alignment(self):
        cases = (
            (">S1\nAAAA\n>S2\nCCC\n", ">S1\nGG\n>S2\nTT\n"),
            (">S1\nAAAA\n>S2\nCCCC\n", ">S1\nGG\n>S2\nT\n"),
        )
        for cpdna_text, mtdna_text in cases:
            with self.subTest(cpdna=cpdna_text, mtdna=mtdna_text):
                with tempfile.TemporaryDirectory() as tmp:
                    root = Path(tmp)
                    cpdna_path = root / "cpDNA.fa"
                    cpdna_path.write_text(cpdna_text)
                    mtdna_path = root / "mtDNA.fa"
                    mtdna_path.write_text(mtdna_text)

                    with self.assertRaisesRegex(
                        ConcatenatedConsensusError,
                        "Inconsistent FASTA sequence lengths",
                    ):
                        concatenate_consensus_alignments(cpdna_path, mtdna_path)

    def test_rejects_empty_fasta_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cpdna_path = root / "cpDNA.fa"
            cpdna_path.write_text("")
            mtdna_path = root / "mtDNA.fa"
            mtdna_path.write_text(">S1\nGG\n")

            with self.assertRaisesRegex(
                ConcatenatedConsensusError,
                "No FASTA records",
            ):
                concatenate_consensus_alignments(cpdna_path, mtdna_path)

    def test_rejects_sequence_before_first_fasta_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cpdna_path = root / "cpDNA.fa"
            cpdna_path.write_text("AAAA\n>S1\nCCCC\n")
            mtdna_path = root / "mtDNA.fa"
            mtdna_path.write_text(">S1\nGG\n")

            with self.assertRaisesRegex(
                ConcatenatedConsensusError,
                "FASTA sequence before header",
            ):
                concatenate_consensus_alignments(cpdna_path, mtdna_path)

    def test_rejects_empty_fasta_sequence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cpdna_path = root / "cpDNA.fa"
            cpdna_path.write_text(">S1\n>S2\n")
            mtdna_path = root / "mtDNA.fa"
            mtdna_path.write_text(">S1\nGG\n>S2\nTT\n")

            with self.assertRaisesRegex(
                ConcatenatedConsensusError,
                "Empty FASTA sequence",
            ):
                concatenate_consensus_alignments(cpdna_path, mtdna_path)


if __name__ == "__main__":
    unittest.main()
