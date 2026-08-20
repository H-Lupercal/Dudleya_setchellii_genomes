import csv
import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.variant_calling import (
    VariantTrack,
    build_bcftools_commands,
    read_variant_samples,
    read_variant_tracks,
    write_variant_call_inputs,
)


class VariantSampleTests(unittest.TestCase):
    def test_read_variant_samples_resolves_safe_bam_paths_and_sample_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample_table = root / "included_samples.tsv"
            bam_dir = root / "bam"
            bam_dir.mkdir()
            sample_table.write_text(
                "sample_id\tbatch\tnaming_profile\tspecies\tpopcode\tpopulation_name\t"
                "du_id\tlp_id\tr1_paths\tr2_paths\tdownstream_cpDNA_use\t"
                "downstream_mtDNA_use\tinclude_reason\n"
                "KEEP/ONE\tbatch\tmain\tD. cymosa\tCY\tPop\tDu-1\tLP_1\t"
                "r1\tr2\tinclude\tinclude\tpasses\n"
                "KEEP_TWO\tbatch\tmain\tD. cymosa\tCY\tPop\tDu-2\tLP_2\t"
                "r1\tr2\tinclude\tinclude\tpasses\n"
            )
            (bam_dir / "KEEP_ONE.organelle.sorted.bam").write_text("")
            (bam_dir / "KEEP_ONE.organelle.sorted.bam.bai").write_text("")
            (bam_dir / "KEEP_TWO.organelle.sorted.bam").write_text("")
            (bam_dir / "KEEP_TWO.organelle.sorted.bam.bai").write_text("")

            samples = read_variant_samples(sample_table, bam_dir, sample_limit=1)

        self.assertEqual(len(samples), 1)
        self.assertEqual(samples[0].sample_id, "KEEP/ONE")
        self.assertEqual(samples[0].safe_sample_id, "KEEP_ONE")
        self.assertEqual(samples[0].bam_path.name, "KEEP_ONE.organelle.sorted.bam")


class VariantTrackTests(unittest.TestCase):
    def test_read_variant_tracks_keeps_only_population_genetic_tracks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            track_table = root / "analysis_tracks.tsv"
            cpdna_bed = root / "cpdna_population_sites.bed"
            mtdna_bed = root / "mtdna_high_confidence_unique_regions.bed"
            cpdna_bed.write_text("chloroplast\t0\t10\tcp\n")
            mtdna_bed.write_text("mitochondria\t0\t5\tmt\n")
            track_table.write_text(
                "track_id\torganelle\tpurpose\tbed_path\tcoordinate_system\tsource\t"
                "step5_use\tnotes\n"
                f"cpdna_population_sites\tcpDNA\tvariant_calling_and_population_genetics\t"
                f"{cpdna_bed}\tBED\tsource\tuse\tnotes\n"
                f"mtdna_high_confidence_unique\tmtDNA\tvariant_calling_and_population_genetics\t"
                f"{mtdna_bed}\tBED\tsource\tuse\tnotes\n"
            )

            tracks = read_variant_tracks(track_table)

        self.assertEqual(set(tracks), {"cpDNA", "mtDNA"})
        self.assertEqual(tracks["cpDNA"].track_id, "cpdna_population_sites")
        self.assertEqual(tracks["mtDNA"].track_id, "mtdna_high_confidence_unique")


class VariantCommandTests(unittest.TestCase):
    def test_write_variant_call_inputs_and_build_commands_use_haploid_track_restricted_calling(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output_dir = root / "08_variant_calling"
            bam_path = root / "sample.bam"
            bam_path.write_text("")
            samples = [
                type(
                    "VariantSample",
                    (),
                    {
                        "sample_id": "SAMPLE_1",
                        "safe_sample_id": "SAMPLE_1",
                        "bam_path": bam_path,
                    },
                )()
            ]
            track = VariantTrack(
                organelle="cpDNA",
                track_id="cpdna_population_sites",
                bed_path=root / "cpdna_population_sites.bed",
                output_prefix="cpDNA.smoke.raw",
            )
            track.bed_path.write_text("chloroplast\t0\t10\tcp\n")

            inputs = write_variant_call_inputs(track, samples, output_dir)
            commands = build_bcftools_commands(
                track=track,
                inputs=inputs,
                reference=Path("reference.fa"),
                min_mapq=20,
                min_baseq=20,
                max_depth=10000,
                threads=4,
            )

            self.assertEqual(inputs.bam_list_path.read_text(), f"{bam_path}\n")
            self.assertEqual(inputs.sample_names_path.read_text(), "SAMPLE_1\n")
            mpileup = commands["mpileup"]
            call = commands["call"]
            self.assertIn("-R", mpileup)
            self.assertIn(str(track.bed_path), mpileup)
            self.assertIn("--ignore-RG", mpileup)
            self.assertIn("--max-depth", mpileup)
            self.assertIn("10000", mpileup)
            self.assertIn("--ploidy", call)
            self.assertIn("1", call)
            self.assertIn("-m", call)
            self.assertIn("-v", call)


if __name__ == "__main__":
    unittest.main()
