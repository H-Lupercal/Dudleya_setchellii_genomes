import csv
import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.analysis_masks import (
    build_cpdna_tracks,
    build_mtdna_tracks,
    generate_analysis_masks,
    interval_to_bed_fields,
)


class AnalysisMaskCoordinateTests(unittest.TestCase):
    def test_interval_to_bed_fields_converts_one_based_inclusive_to_zero_based_bed(self):
        fields = interval_to_bed_fields(
            record="chloroplast",
            start_1based=82091,
            end_1based=107826,
            name="cpdna_IR_copy_1",
        )

        self.assertEqual(fields, ["chloroplast", "82090", "107826", "cpdna_IR_copy_1"])


class ChloroplastMaskTests(unittest.TestCase):
    def test_cpdna_tracks_keep_one_ir_copy_for_population_sites(self):
        with tempfile.TemporaryDirectory() as tmp:
            repeat_path = Path(tmp) / "cpdna_self_repeat_intervals.tsv"
            repeat_path.write_text(
                "rank\tlength_bp\tidentity_percent\tquery_start\tquery_end\t"
                "match_start\tmatch_end\torientation\n"
                "1\t4\t99.9\t15\t18\t5\t8\treverse\n"
                "2\t4\t99.9\t5\t8\t15\t18\treverse\n"
            )

            tracks = build_cpdna_tracks(repeat_path, reference_length=20)

        ir_regions = [(row.start_1based, row.end_1based, row.name) for row in tracks.ir_regions]
        population_sites = [
            (row.start_1based, row.end_1based, row.name)
            for row in tracks.population_regions
        ]
        duplicate_mask = [
            (row.start_1based, row.end_1based, row.name)
            for row in tracks.duplicate_ir_mask
        ]

        self.assertEqual(
            ir_regions,
            [(5, 8, "cpdna_IR_copy_1"), (15, 18, "cpdna_IR_copy_2")],
        )
        self.assertEqual(duplicate_mask, [(15, 18, "cpdna_duplicate_IR_copy_mask")])
        self.assertEqual(
            population_sites,
            [
                (1, 14, "cpdna_population_single_IR_region_1"),
                (19, 20, "cpdna_population_single_IR_region_2"),
            ],
        )


class MitochondrialMaskTests(unittest.TestCase):
    def test_mtdna_tracks_split_permissive_coverage_and_high_confidence_unique_regions(self):
        with tempfile.TemporaryDirectory() as tmp:
            interval_path = Path(tmp) / "mtdna_high_mapq_consensus_intervals.tsv"
            interval_path.write_text(
                "threshold_usable_samples\trank\tstart\tend\tlength_bp\n"
                "12\t1\t4\t8\t5\n"
                "12\t2\t20\t25\t6\n"
                "10\t1\t3\t9\t7\n"
            )

            tracks = build_mtdna_tracks(
                interval_path,
                reference_length=30,
                high_confidence_threshold=12,
            )

        permissive = [
            (row.start_1based, row.end_1based, row.name)
            for row in tracks.permissive_coverage_regions
        ]
        unique = [
            (row.start_1based, row.end_1based, row.name)
            for row in tracks.high_confidence_unique_regions
        ]

        self.assertEqual(permissive, [(1, 30, "mtdna_permissive_coverage_full_reference")])
        self.assertEqual(
            unique,
            [
                (4, 8, "mtdna_high_confidence_unique_region_1"),
                (20, 25, "mtdna_high_confidence_unique_region_2"),
            ],
        )


class AnalysisMaskOutputTests(unittest.TestCase):
    def test_generate_analysis_masks_writes_bed_tracks_and_track_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cpdna_repeats = root / "cpdna_self_repeat_intervals.tsv"
            mtdna_intervals = root / "mtdna_high_mapq_consensus_intervals.tsv"
            output_dir = root / "05_analysis_masks"
            cpdna_repeats.write_text(
                "rank\tlength_bp\tidentity_percent\tquery_start\tquery_end\t"
                "match_start\tmatch_end\torientation\n"
                "1\t4\t99.9\t15\t18\t5\t8\treverse\n"
                "2\t4\t99.9\t5\t8\t15\t18\treverse\n"
            )
            mtdna_intervals.write_text(
                "threshold_usable_samples\trank\tstart\tend\tlength_bp\n"
                "12\t1\t4\t8\t5\n"
                "12\t2\t20\t25\t6\n"
            )

            generate_analysis_masks(
                output_dir=output_dir,
                cpdna_repeat_intervals_path=cpdna_repeats,
                mtdna_high_mapq_intervals_path=mtdna_intervals,
                cpdna_length=20,
                mtdna_length=30,
                mtdna_high_confidence_threshold=12,
            )

            with (output_dir / "analysis_tracks.tsv").open(newline="") as handle:
                tracks = list(csv.DictReader(handle, delimiter="\t"))

            self.assertTrue((output_dir / "cpdna_ir_regions.bed").exists())
            self.assertTrue((output_dir / "cpdna_population_sites.bed").exists())
            self.assertTrue(
                (output_dir / "mtdna_permissive_coverage_regions.bed").exists()
            )
            self.assertTrue(
                (output_dir / "mtdna_high_confidence_unique_regions.bed").exists()
            )
            self.assertTrue((output_dir / "mask_summary.md").exists())
            self.assertEqual(
                [row["track_id"] for row in tracks],
                [
                    "cpdna_full_coverage",
                    "cpdna_ir_regions",
                    "cpdna_duplicate_ir_mask",
                    "cpdna_population_sites",
                    "mtdna_permissive_coverage",
                    "mtdna_high_confidence_unique",
                ],
            )


if __name__ == "__main__":
    unittest.main()
