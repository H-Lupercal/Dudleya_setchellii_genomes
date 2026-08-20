import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.admixture_analysis import (
    AdmixtureInput,
    build_admixture_command,
    parse_cv_error,
    read_admixture_inputs,
    summarize_replicate_stability,
    write_pseudo_diploid_ped_map,
    write_admixture_outputs,
)


class AdmixtureInputTests(unittest.TestCase):
    def test_read_admixture_inputs_uses_snp_alignment_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            snp_dir = Path(tmp)
            fasta_path = snp_dir / "cpDNA.primary.snp_alignment.fa"
            site_path = snp_dir / "cpDNA.primary.snp_sites.tsv"
            fasta_path.write_text(">DU-1\nAC\n>DU-2\nGT\n")
            site_path.write_text("organelle\tposition\ncpDNA\t1\ncpDNA\t2\n")
            (snp_dir / "primary.snp_alignment_summary.tsv").write_text(
                "organelle\ttrack_id\tsample_count\tfiltered_records\talignment_sites\t"
                "missing_bases\tfiltered_vcf_path\talignment_fasta_path\tsite_table_path\n"
                f"cpDNA\tcpdna_population_sites\t2\t2\t2\t0\tcp.vcf.gz\t"
                f"{fasta_path}\t{site_path}\n"
            )

            inputs = read_admixture_inputs(snp_dir, run_label="primary")

        self.assertEqual(len(inputs), 1)
        self.assertEqual(inputs[0].organelle, "cpDNA")
        self.assertEqual(inputs[0].alignment_fasta_path, fasta_path)

    def test_write_pseudo_diploid_ped_map_duplicates_haploid_calls(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            fasta_path = output_dir / "tiny.fa"
            fasta_path.write_text(">DU-1\nAN\n>DU-2\nTG\n")
            admixture_input = AdmixtureInput(
                organelle="mtDNA",
                track_id="mtdna_high_confidence_unique",
                sample_count=2,
                alignment_sites=2,
                missing_bases=1,
                alignment_fasta_path=fasta_path,
                site_table_path=Path("sites.tsv"),
            )

            ped_path, map_path, included_sample_ids, excluded_sample_ids = (
                write_pseudo_diploid_ped_map(
                    admixture_input,
                    output_dir,
                    run_label="primary",
                )
            )

            ped = ped_path.read_text()
            map_text = map_path.read_text()

        self.assertEqual(included_sample_ids, ["DU-1", "DU-2"])
        self.assertEqual(excluded_sample_ids, [])
        self.assertIn("DU-1 DU-1 0 0 0 -9 A A 0 0", ped)
        self.assertIn("DU-2 DU-2 0 0 0 -9 T T G G", ped)
        self.assertIn("mtDNA_snp_1", map_text)

    def test_write_pseudo_diploid_ped_map_excludes_all_missing_samples(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            fasta_path = output_dir / "tiny.fa"
            fasta_path.write_text(">DU-1\nNN\n>DU-2\nTG\n")
            admixture_input = AdmixtureInput(
                organelle="mtDNA",
                track_id="mtdna_high_confidence_unique",
                sample_count=2,
                alignment_sites=2,
                missing_bases=2,
                alignment_fasta_path=fasta_path,
                site_table_path=Path("sites.tsv"),
            )

            ped_path, _, included_sample_ids, excluded_sample_ids = write_pseudo_diploid_ped_map(
                admixture_input,
                output_dir,
                run_label="primary",
            )

            ped = ped_path.read_text()
            excluded = (
                output_dir
                / "mtDNA.primary.pseudo_diploid.excluded_samples.tsv"
            ).read_text()

        self.assertEqual(included_sample_ids, ["DU-2"])
        self.assertEqual(excluded_sample_ids, ["DU-1"])
        self.assertNotIn("DU-1 DU-1", ped)
        self.assertIn("DU-2 DU-2 0 0 0 -9 T T G G", ped)
        self.assertIn("DU-1\tmtDNA\tall_snp_genotypes_missing", excluded)


class AdmixtureCommandTests(unittest.TestCase):
    def test_build_admixture_command_uses_cv_seed_threads_and_k(self):
        command = build_admixture_command("admixture", Path("cpDNA.ped"), k=3, seed=42, threads=4)

        self.assertEqual(command, ["admixture", "--cv", "--seed=42", "-j4", "cpDNA.ped", "3"])

    def test_parse_cv_error_reads_admixture_stdout(self):
        log = "Random text\nCV error (K=4): 0.12345\n"

        self.assertAlmostEqual(parse_cv_error(log), 0.12345)

    def test_summarize_replicate_stability_selects_best_mean_cv_and_counts_replicates(self):
        rows = [
            {"organelle": "cpDNA", "k": "2", "replicate": "1", "cv_error": "0.30"},
            {"organelle": "cpDNA", "k": "2", "replicate": "2", "cv_error": "0.32"},
            {"organelle": "cpDNA", "k": "3", "replicate": "1", "cv_error": "0.20"},
            {"organelle": "cpDNA", "k": "3", "replicate": "2", "cv_error": "0.22"},
        ]

        summary = summarize_replicate_stability(rows)

        self.assertEqual(summary[0]["k"], "2")
        self.assertEqual(summary[0]["replicate_count"], "2")
        self.assertEqual(summary[1]["is_best_mean_k"], "yes")
        self.assertEqual(summary[1]["mean_cv_error"], "0.21000000")


class AdmixtureOutputTests(unittest.TestCase):
    def test_write_admixture_outputs_records_summary_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            rows = [
                {
                    "organelle": "cpDNA",
                    "track_id": "cpdna_population_sites",
                    "k": "2",
                    "replicate": "1",
                    "cv_error": "0.1",
                    "is_best_k": "yes",
                    "is_best_mean_k": "yes",
                    "mean_cv_error": "0.1",
                    "sd_cv_error": "0.0",
                    "replicate_count": "1",
                    "excluded_sample_count": "0",
                    "q_path": "cpDNA.2.Q",
                    "p_path": "cpDNA.2.P",
                    "log_path": "cpDNA.K2.log",
                    "best_q_table_path": "cpDNA.bestK2.q.tsv",
                    "structure_png_path": "cpDNA.bestK2.structure.png",
                    "structure_pdf_path": "cpDNA.bestK2.structure.pdf",
                    "structure_svg_path": "cpDNA.bestK2.structure.svg",
                    "cv_plot_path": "cpDNA.admixture_cv.png",
                    "plink_command": "plink --file cpDNA --make-bed",
                    "command": "admixture --cv cpDNA.ped 2",
                }
            ]

            write_admixture_outputs(output_dir, rows, run_label="primary")

            summary = (output_dir / "primary.admixture_summary.tsv").read_text()
            report = (output_dir / "primary.admixture_report.md").read_text()

        self.assertIn("cv_error", summary)
        self.assertIn("# Admixture-Style Clustering", report)
        self.assertIn("pseudo-diploid homozygous", report)


if __name__ == "__main__":
    unittest.main()
