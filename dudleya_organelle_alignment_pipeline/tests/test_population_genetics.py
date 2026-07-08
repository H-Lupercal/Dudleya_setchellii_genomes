import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.population_genetics import (
    compute_haplotype_diversity,
    compute_pairwise_fst,
    read_population_inputs,
    run_one_population_summary,
    write_population_genetics_outputs,
)


class PopulationGeneticsInputTests(unittest.TestCase):
    def test_read_population_inputs_uses_snp_alignment_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            snp_dir = Path(tmp)
            fasta_path = snp_dir / "cpDNA.primary.snp_alignment.fa"
            site_path = snp_dir / "cpDNA.primary.snp_sites.tsv"
            fasta_path.write_text(">S1\nAA\n>S2\nTT\n")
            site_path.write_text("organelle\tposition\ncpDNA\t1\n")
            (snp_dir / "primary.snp_alignment_summary.tsv").write_text(
                "organelle\ttrack_id\tsample_count\tfiltered_records\talignment_sites\t"
                "missing_bases\tfiltered_vcf_path\talignment_fasta_path\tsite_table_path\n"
                f"cpDNA\tcpdna_population_sites\t2\t2\t2\t0\tcp.vcf.gz\t"
                f"{fasta_path}\t{site_path}\n"
            )

            inputs = read_population_inputs(snp_dir, run_label="primary")

        self.assertEqual(inputs[0].organelle, "cpDNA")
        self.assertEqual(inputs[0].alignment_sites, 2)


class PopulationMetricTests(unittest.TestCase):
    def test_compute_haplotype_diversity_counts_unique_sequences(self):
        diversity = compute_haplotype_diversity(["AA", "AA", "AT"])

        self.assertAlmostEqual(diversity, 2 / 3)

    def test_compute_pairwise_fst_returns_high_value_for_fixed_differences(self):
        fst, informative_sites = compute_pairwise_fst(["AA", "AA"], ["TT", "TT"])

        self.assertAlmostEqual(fst, 1.0)
        self.assertEqual(informative_sites, 2)


class PopulationOutputTests(unittest.TestCase):
    def test_run_one_population_summary_writes_pairwise_and_population_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp) / "out"
            metadata_path = Path(tmp) / "included_samples.tsv"
            fasta_path = Path(tmp) / "mtDNA.primary.snp_alignment.fa"
            metadata_path.write_text(
                "sample_id\tspecies\tpopcode\tpopulation_name\tnaming_profile\n"
                "S1\tCY\tRED\tRed Pop\tmain_standard\n"
                "S2\tCY\tRED\tRed Pop\tmain_standard\n"
                "S3\tCY\tBLU\tBlue Pop\tmain_standard\n"
                "S4\tCY\tBLU\tBlue Pop\tmain_standard\n"
            )
            fasta_path.write_text(">S1\nAA\n>S2\nAT\n>S3\nTT\n>S4\nTT\n")

            result = run_one_population_summary(
                organelle="mtDNA",
                track_id="mtdna_high_confidence_unique",
                alignment_sites=2,
                alignment_fasta_path=fasta_path,
                metadata_path=metadata_path,
                output_dir=output_dir,
                run_label="primary",
            )
            write_population_genetics_outputs(output_dir, [result], run_label="primary")

            pairwise = result.pairwise_fst_path.read_text()
            population = result.population_summary_path.read_text()
            report = (output_dir / "primary.population_genetics_report.md").read_text()

        self.assertIn("fst", pairwise)
        self.assertIn("haplotype_diversity", population)
        self.assertIn("# Step 17 Population Genetics", report)


if __name__ == "__main__":
    unittest.main()
