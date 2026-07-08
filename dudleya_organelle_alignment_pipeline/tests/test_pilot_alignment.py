import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.pilot_alignment import (
    AlignmentError,
    OrganelleMetrics,
    build_depth_command,
    build_sample_summary,
    count_fastq_records,
    parse_depth_file,
    parse_idxstats_file,
    read_alignment_samples,
    safe_sample_name,
)


class PilotAlignmentMetricsTests(unittest.TestCase):
    def test_depth_parser_counts_missing_positions_as_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            depth_path = Path(tmp) / "depth.tsv"
            depth_path.write_text(
                "chloroplast\t1\t3\n"
                "chloroplast\t2\t7\n"
                "mitochondria\t1\t12\n"
            )

            metrics = parse_depth_file(
                depth_path,
                {"chloroplast": 4, "mitochondria": 2},
            )

        self.assertEqual(metrics["chloroplast"].total_depth, 10)
        self.assertEqual(metrics["chloroplast"].bases_ge_1x, 2)
        self.assertEqual(metrics["chloroplast"].bases_ge_5x, 1)
        self.assertEqual(metrics["chloroplast"].bases_ge_10x, 0)
        self.assertEqual(metrics["chloroplast"].mean_depth, 2.5)
        self.assertEqual(metrics["chloroplast"].breadth_ge_1x, 0.5)
        self.assertEqual(metrics["mitochondria"].total_depth, 12)
        self.assertEqual(metrics["mitochondria"].bases_ge_10x, 1)
        self.assertEqual(metrics["mitochondria"].breadth_ge_10x, 0.5)

    def test_idxstats_parser_keeps_reference_mapped_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            idxstats_path = Path(tmp) / "sample.idxstats.tsv"
            idxstats_path.write_text(
                "chloroplast\t150274\t42\t3\n"
                "mitochondria\t243359\t7\t1\n"
                "*\t0\t0\t9\n"
            )

            counts = parse_idxstats_file(idxstats_path)

        self.assertEqual(counts, {"chloroplast": 42, "mitochondria": 7})

    def test_sample_summary_flags_unbalanced_cross_organelle_signal(self):
        summary = build_sample_summary(
            sample_id="S1",
            row={"batch": "batch", "species": "D. setchellii", "popcode": "BAI"},
            mapped_counts={"chloroplast": 1000, "mitochondria": 2},
            depth_metrics={
                "chloroplast": OrganelleMetrics("chloroplast", 100, 1000, 90, 80, 70),
                "mitochondria": OrganelleMetrics("mitochondria", 100, 10, 2, 1, 0),
            },
        )

        self.assertEqual(summary["total_organelle_mapped_reads"], "1002")
        self.assertEqual(summary["chloroplast_fraction_of_organelle_mapped"], "0.998004")
        self.assertIn("low_mitochondria_mapped_reads", summary["qc_notes"])

    def test_safe_sample_name_removes_path_unfriendly_characters(self):
        self.assertEqual(safe_sample_name("CY_ALA_LP_298-Du-767"), "CY_ALA_LP_298-Du-767")
        self.assertEqual(safe_sample_name("odd sample/1"), "odd_sample_1")

    def test_count_fastq_records_counts_four_line_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            fastq_path = Path(tmp) / "reads.fastq"
            fastq_path.write_text("@r1\nAC\n+\nII\n@r2\nTG\n+\nII\n")

            count = count_fastq_records(fastq_path)

        self.assertEqual(count, 2)

    def test_depth_command_uses_samtools_quality_flags_correctly(self):
        command = build_depth_command(
            Path("sample.bam"),
            min_mapq=7,
            min_baseq=19,
        )

        self.assertEqual(
            command,
            [
                "samtools",
                "depth",
                "-aa",
                "-q",
                "19",
                "-Q",
                "7",
                "sample.bam",
            ],
        )


class PilotAlignmentSampleTableTests(unittest.TestCase):
    def test_read_alignment_samples_rejects_missing_mates(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "pilot.tsv"
            table.write_text(
                "sample_id\tr1_paths\tr2_paths\tanalysis_status\tpair_status\n"
                "ok\treads_R1.fastq.gz\treads_R2.fastq.gz\tinclude_primary_paired_end\tcomplete\n"
                "missing\t\treads_R2.fastq.gz\texclude_missing_mate\tmissing_R1\n"
            )

            samples = read_alignment_samples(table)

        self.assertEqual([sample.sample_id for sample in samples], ["ok"])

    def test_read_alignment_samples_errors_on_multi_file_rows_for_now(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = Path(tmp) / "pilot.tsv"
            table.write_text(
                "sample_id\tr1_paths\tr2_paths\tanalysis_status\tpair_status\n"
                "multi\tr1a.fastq.gz;r1b.fastq.gz\tr2a.fastq.gz;r2b.fastq.gz\t"
                "include_primary_paired_end\tcomplete\n"
            )

            with self.assertRaises(AlignmentError):
                read_alignment_samples(table)


if __name__ == "__main__":
    unittest.main()
