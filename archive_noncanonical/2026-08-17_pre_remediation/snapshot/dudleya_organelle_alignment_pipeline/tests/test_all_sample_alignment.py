import csv
import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.all_sample_alignment import (
    build_track_summary_rows,
    parse_track_depth_file,
    read_track_regions,
)


class AllSampleTrackParsingTests(unittest.TestCase):
    def test_read_track_regions_loads_bed_coordinates_from_track_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bed = root / "cpdna_population_sites.bed"
            bed.write_text("chloroplast\t0\t3\tcp_segment_1\n")
            tracks = root / "analysis_tracks.tsv"
            tracks.write_text(
                "track_id\torganelle\tpurpose\tbed_path\tcoordinate_system\t"
                "source\tstep5_use\tnotes\n"
                f"cpdna_population_sites\tcpDNA\tvariant_calling_and_population_genetics\t{bed}\t"
                "BED 0-based half-open\tsource\tuse\tnotes\n"
            )

            regions = read_track_regions(tracks)

        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].track_id, "cpdna_population_sites")
        self.assertEqual(regions[0].record, "chloroplast")
        self.assertEqual(regions[0].start_1based, 1)
        self.assertEqual(regions[0].end_1based, 3)
        self.assertEqual(regions[0].length_bp, 3)


class AllSampleTrackDepthTests(unittest.TestCase):
    def test_parse_track_depth_file_counts_only_sites_inside_track_regions(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cp_bed = root / "cp.bed"
            mt_bed = root / "mt.bed"
            cp_bed.write_text("chloroplast\t0\t3\tcp_region\n")
            mt_bed.write_text("mitochondria\t1\t4\tmt_region\n")
            tracks = root / "analysis_tracks.tsv"
            tracks.write_text(
                "track_id\torganelle\tpurpose\tbed_path\tcoordinate_system\t"
                "source\tstep5_use\tnotes\n"
                f"cpdna_population_sites\tcpDNA\tvariant_calling_and_population_genetics\t{cp_bed}\t"
                "BED 0-based half-open\tsource\tuse\tnotes\n"
                f"mtdna_high_confidence_unique\tmtDNA\tvariant_calling_and_population_genetics\t{mt_bed}\t"
                "BED 0-based half-open\tsource\tuse\tnotes\n"
            )
            depth = root / "sample.depth.tsv"
            depth.write_text(
                "chloroplast\t1\t3\n"
                "chloroplast\t2\t0\n"
                "chloroplast\t3\t6\n"
                "chloroplast\t4\t20\n"
                "mitochondria\t1\t30\n"
                "mitochondria\t3\t7\n"
            )

            metrics = parse_track_depth_file(depth, read_track_regions(tracks))

        cp_metrics = metrics["cpdna_population_sites"]
        mt_metrics = metrics["mtdna_high_confidence_unique"]
        self.assertEqual(cp_metrics.region_bp, 3)
        self.assertEqual(cp_metrics.total_depth, 9)
        self.assertEqual(cp_metrics.bases_ge_1x, 2)
        self.assertEqual(cp_metrics.bases_ge_5x, 1)
        self.assertEqual(cp_metrics.bases_ge_10x, 0)
        self.assertEqual(mt_metrics.region_bp, 3)
        self.assertEqual(mt_metrics.total_depth, 7)
        self.assertEqual(mt_metrics.bases_ge_1x, 1)
        self.assertEqual(mt_metrics.bases_ge_5x, 1)
        self.assertEqual(mt_metrics.bases_ge_10x, 0)

    def test_build_track_summary_rows_includes_track_purpose_and_breadth(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bed = root / "cp.bed"
            bed.write_text("chloroplast\t0\t3\tcp_region\n")
            tracks = root / "analysis_tracks.tsv"
            tracks.write_text(
                "track_id\torganelle\tpurpose\tbed_path\tcoordinate_system\t"
                "source\tstep5_use\tnotes\n"
                f"cpdna_population_sites\tcpDNA\tvariant_calling_and_population_genetics\t{bed}\t"
                "BED 0-based half-open\tsource\tuse\tnotes\n"
            )
            depth = root / "sample.depth.tsv"
            depth.write_text("chloroplast\t1\t3\nchloroplast\t3\t6\n")
            regions = read_track_regions(tracks)
            metrics = parse_track_depth_file(depth, regions)

            rows = build_track_summary_rows(
                sample_id="S1",
                row={"batch": "batch", "species": "D. setchellii", "popcode": "BAI"},
                track_regions=regions,
                track_metrics=metrics,
                depth_path=depth,
            )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["track_id"], "cpdna_population_sites")
        self.assertEqual(rows[0]["purpose"], "variant_calling_and_population_genetics")
        self.assertEqual(rows[0]["region_bp"], "3")
        self.assertEqual(rows[0]["breadth_ge_1x"], "0.666667")
        self.assertEqual(rows[0]["depth_path"], depth.as_posix())


if __name__ == "__main__":
    unittest.main()
