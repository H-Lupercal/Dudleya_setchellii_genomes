import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.pca_analysis import (
    PcaInput,
    build_haploid_snp_matrix,
    read_pca_inputs,
    run_one_pca,
    write_pca_outputs,
)


class PcaInputTests(unittest.TestCase):
    def test_read_pca_inputs_uses_snp_alignment_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            snp_dir = Path(tmp)
            fasta_path = snp_dir / "cpDNA.primary.snp_alignment.fa"
            site_path = snp_dir / "cpDNA.primary.snp_sites.tsv"
            fasta_path.write_text(">DU-1\nACGT\n>DU-2\nTCGA\n")
            site_path.write_text("site\n1\n")
            (snp_dir / "primary.snp_alignment_summary.tsv").write_text(
                "organelle\ttrack_id\tsample_count\tfiltered_records\talignment_sites\t"
                "missing_bases\tfiltered_vcf_path\talignment_fasta_path\tsite_table_path\n"
                f"cpDNA\tcpdna_population_sites\t2\t4\t4\t0\tcp.vcf.gz\t"
                f"{fasta_path}\t{site_path}\n"
            )

            inputs = read_pca_inputs(snp_dir, run_label="primary")

        self.assertEqual(len(inputs), 1)
        self.assertEqual(inputs[0].organelle, "cpDNA")
        self.assertEqual(inputs[0].alignment_sites, 4)
        self.assertEqual(inputs[0].alignment_fasta_path, fasta_path)


class PcaMatrixTests(unittest.TestCase):
    def test_build_haploid_snp_matrix_encodes_bases_and_imputes_missing_site_mean(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "tiny.fa"
            fasta.write_text(">DU-1\nAAN\n>DU-2\nATG\n>DU-3\nTTT\n")

            matrix, sample_ids, retained_sites = build_haploid_snp_matrix(fasta)

        self.assertEqual(sample_ids, ["DU-1", "DU-2", "DU-3"])
        self.assertEqual(matrix.shape, (3, 3))
        self.assertEqual(retained_sites, 3)
        self.assertAlmostEqual(float(matrix[0, 2]), 0.5)


class PcaOutputTests(unittest.TestCase):
    def test_run_one_pca_writes_coordinates_variance_and_figures(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp) / "out"
            metadata_path = Path(tmp) / "included_samples.tsv"
            fasta_path = Path(tmp) / "mtDNA.primary.snp_alignment.fa"
            metadata_path.write_text(
                "sample_id\tspecies\tpopcode\tpopulation_name\tnaming_profile\n"
                "DU-1\tCY\tRED\tRed Pop\tmain_standard\n"
                "DU-2\tCY\tRED\tRed Pop\tmain_standard\n"
                "DU-3\tAB\tMAD\tMad Pop\tmain_standard\n"
                "DU-4\tAB\tMAD\tMad Pop\tmain_standard\n"
            )
            fasta_path.write_text(
                ">DU-1\nAAAA\n"
                ">DU-2\nAAAT\n"
                ">DU-3\nTTTA\n"
                ">DU-4\nTTTT\n"
            )
            pca_input = PcaInput(
                organelle="mtDNA",
                track_id="mtdna_high_confidence_unique",
                sample_count=4,
                alignment_sites=4,
                missing_bases=0,
                alignment_fasta_path=fasta_path,
                site_table_path=Path("sites.tsv"),
            )

            result = run_one_pca(pca_input, metadata_path, output_dir, run_label="primary")

            write_pca_outputs(output_dir, [result], run_label="primary")

            summary = (output_dir / "primary.pca_summary.tsv").read_text()
            report = (output_dir / "primary.pca_report.md").read_text()

            self.assertTrue(result.coordinates_path.exists())
            self.assertTrue(result.variance_path.exists())
            self.assertTrue(result.png_path.exists())
            self.assertTrue(result.pdf_path.exists())
            self.assertTrue(result.svg_path.exists())
            self.assertIn("mtDNA", summary)
            self.assertIn("# PCA Visualization", report)


if __name__ == "__main__":
    unittest.main()
