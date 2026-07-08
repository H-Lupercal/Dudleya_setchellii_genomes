import csv
import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.downstream_sample_set import (
    build_downstream_sample_set,
    write_downstream_sample_set_outputs,
)


def manifest_row(sample_id: str, status: str = "include_primary_paired_end") -> dict[str, str]:
    return {
        "sample_id": sample_id,
        "batch": "batch",
        "naming_profile": "main_standard",
        "popcode": sample_id.split("_LP_")[0] if "_LP_" in sample_id else "",
        "species": "D. cymosa",
        "population_name": "Population",
        "du_id": "Du-1",
        "lp_id": "LP_1",
        "sequencing_samples": "S1",
        "lanes": "L005",
        "r1_paths": f"{sample_id}_R1.fastq.gz",
        "r2_paths": f"{sample_id}_R2.fastq.gz",
        "r1_count": "1",
        "r2_count": "1",
        "pair_status": "complete" if status == "include_primary_paired_end" else "missing_R2",
        "metadata_status": "resolved",
        "analysis_status": status,
        "analysis_note": "Use in the primary paired-end cpDNA/mtDNA alignment workflow.",
    }


class DownstreamSampleSetTests(unittest.TestCase):
    def test_build_downstream_sample_set_excludes_qc_ignored_and_missing_mates(self):
        analysis_rows = [
            manifest_row("KEEP_LP_001_Du-1"),
            manifest_row("DROP_QC_LP_002_Du-2"),
            manifest_row("KEEP_LP_003_Du-3"),
        ]
        upstream_excluded_rows = [
            manifest_row("MISSING_LP_004_Du-4", status="exclude_missing_mate")
        ]
        qc_decision_rows = [
            {
                "sample_id": "DROP_QC_LP_002_Du-2",
                "downstream_cpDNA_use": "exclude",
                "downstream_mtDNA_use": "exclude",
                "ignored_downstream": "yes",
                "reason": "tiny FASTQ input",
                "evidence": "input_read_records=42",
            }
        ]

        included, excluded = build_downstream_sample_set(
            analysis_rows=analysis_rows,
            upstream_excluded_rows=upstream_excluded_rows,
            qc_decision_rows=qc_decision_rows,
        )

        self.assertEqual(
            [row["sample_id"] for row in included],
            ["KEEP_LP_001_Du-1", "KEEP_LP_003_Du-3"],
        )
        self.assertEqual(
            [row["sample_id"] for row in excluded],
            ["DROP_QC_LP_002_Du-2", "MISSING_LP_004_Du-4"],
        )
        self.assertEqual(excluded[0]["exclusion_stage"], "step5_downstream_qc")
        self.assertEqual(excluded[0]["exclusion_reason"], "tiny FASTQ input")
        self.assertEqual(excluded[1]["exclusion_stage"], "step0_manifest")
        self.assertEqual(excluded[1]["exclusion_reason"], "missing_R2")

    def test_write_downstream_sample_set_outputs_writes_three_expected_files(self):
        included = [
            {
                "sample_id": "KEEP_LP_001_Du-1",
                "batch": "batch",
                "naming_profile": "main_standard",
                "species": "D. cymosa",
                "popcode": "KEEP",
                "population_name": "Population",
                "du_id": "Du-1",
                "lp_id": "LP_1",
                "r1_paths": "r1.fastq.gz",
                "r2_paths": "r2.fastq.gz",
                "downstream_cpDNA_use": "include",
                "downstream_mtDNA_use": "include",
                "include_reason": "passes Step 5 downstream QC",
            }
        ]
        excluded = [
            {
                "sample_id": "DROP_QC_LP_002_Du-2",
                "batch": "batch",
                "naming_profile": "main_standard",
                "species": "D. cymosa",
                "popcode": "DROP_QC",
                "population_name": "Population",
                "du_id": "Du-2",
                "lp_id": "LP_2",
                "exclusion_stage": "step5_downstream_qc",
                "downstream_cpDNA_use": "exclude",
                "downstream_mtDNA_use": "exclude",
                "exclusion_reason": "tiny FASTQ input",
                "evidence": "input_read_records=42",
            }
        ]

        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            write_downstream_sample_set_outputs(included, excluded, output_dir)
            with (output_dir / "included_samples.tsv").open(newline="") as handle:
                included_rows = list(csv.DictReader(handle, delimiter="\t"))
            with (output_dir / "excluded_samples.tsv").open(newline="") as handle:
                excluded_rows = list(csv.DictReader(handle, delimiter="\t"))
            report_text = (output_dir / "downstream_sample_set_report.md").read_text()

        self.assertEqual(included_rows[0]["sample_id"], "KEEP_LP_001_Du-1")
        self.assertEqual(excluded_rows[0]["sample_id"], "DROP_QC_LP_002_Du-2")
        self.assertIn("Included samples: 1", report_text)
        self.assertIn("Excluded samples: 1", report_text)


if __name__ == "__main__":
    unittest.main()
