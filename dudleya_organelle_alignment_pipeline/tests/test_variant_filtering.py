import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.variant_filtering import (
    FilterInput,
    FilterResult,
    build_filter_command,
    write_filtering_outputs,
)


class VariantFilteringCommandTests(unittest.TestCase):
    def test_build_filter_command_keeps_biallelic_snps_and_missingness_threshold(self):
        raw_vcf = Path("cpDNA.primary.raw.vcf.gz")
        filtered_vcf = Path("cpDNA.primary.filtered.vcf.gz")

        command = build_filter_command(
            raw_vcf=raw_vcf,
            filtered_vcf=filtered_vcf,
            max_missing_fraction=0.2,
            min_minor_allele_count=2,
            threads=4,
        )

        self.assertEqual(command[:2], ["bcftools", "view"])
        self.assertIn("-m2", command)
        self.assertIn("-M2", command)
        self.assertIn("-v", command)
        self.assertIn("snps", command)
        self.assertIn("--min-ac", command)
        self.assertIn("2:minor", command)
        self.assertIn("-i", command)
        self.assertIn("F_MISSING<=0.2", command)
        self.assertIn("-Oz", command)
        self.assertIn(str(filtered_vcf), command)
        self.assertIn(str(raw_vcf), command)


class VariantFilteringOutputTests(unittest.TestCase):
    def test_write_filtering_outputs_records_summary_commands_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            result = FilterResult(
                organelle="cpDNA",
                track_id="cpdna_population_sites",
                sample_count=275,
                raw_records=2475,
                filtered_records=2000,
                raw_vcf_path=Path("cpDNA.primary.raw.vcf.gz"),
                filtered_vcf_path=Path("cpDNA.primary.filtered.vcf.gz"),
                filtered_index_path=Path("cpDNA.primary.filtered.vcf.gz.tbi"),
                log_path=Path("cpDNA.primary.filtered.bcftools.log"),
            )

            write_filtering_outputs(
                output_dir=output_dir,
                results=[result],
                command_rows=[
                    {
                        "organelle": "cpDNA",
                        "track_id": "cpdna_population_sites",
                        "step": "filter",
                        "command": "bcftools view ...",
                    }
                ],
                run_label="primary",
                max_missing_fraction=0.2,
                min_minor_allele_count=2,
            )

            summary = (output_dir / "primary.variant_filtering_summary.tsv").read_text()
            report = (output_dir / "primary.variant_filtering_report.md").read_text()
            commands = (output_dir / "primary.filtering.commands.tsv").read_text()

        self.assertIn("filtered_records", summary)
        self.assertIn("2000", summary)
        self.assertIn("# Step 8 Variant Filtering", report)
        self.assertIn("Maximum missing genotype fraction: 0.2", report)
        self.assertIn("Minimum minor allele count: 2", report)
        self.assertIn("bcftools view", commands)


class FilterInputTests(unittest.TestCase):
    def test_filter_input_derives_primary_filtered_paths_from_raw_paths(self):
        filter_input = FilterInput(
            organelle="mtDNA",
            track_id="mtdna_high_confidence_unique",
            sample_count=275,
            raw_records=190,
            raw_vcf_path=Path("mtDNA.primary.raw.vcf.gz"),
            raw_vcf_index_path=Path("mtDNA.primary.raw.vcf.gz.tbi"),
        )

        self.assertEqual(filter_input.filtered_vcf_path.name, "mtDNA.primary.filtered.vcf.gz")
        self.assertEqual(
            filter_input.filtered_index_path.name,
            "mtDNA.primary.filtered.vcf.gz.tbi",
        )
        self.assertEqual(
            filter_input.log_path.name,
            "mtDNA.primary.filtered.bcftools.log",
        )


if __name__ == "__main__":
    unittest.main()
