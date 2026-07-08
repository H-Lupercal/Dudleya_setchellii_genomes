import gzip
import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.snp_alignment import (
    SnpAlignmentInput,
    build_snp_alignment,
    write_alignment_outputs,
)


class SnpAlignmentBuildTests(unittest.TestCase):
    def test_build_snp_alignment_encodes_haploid_genotypes_as_bases(self):
        with tempfile.TemporaryDirectory() as tmp:
            vcf_path = Path(tmp) / "tiny.vcf.gz"
            with gzip.open(vcf_path, "wt") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                    "SAMPLE_REF\tSAMPLE_ALT\tSAMPLE_MISSING\n"
                )
                handle.write(
                    "chloroplast\t10\t.\tA\tG\t.\tPASS\t.\tGT:DP\t"
                    "0:12\t1:9\t.:0\n"
                )
                handle.write(
                    "chloroplast\t20\t.\tT\tC\t.\tPASS\t.\tGT:DP\t"
                    "1:8\t0:11\t./.:0\n"
                )

            alignment = build_snp_alignment(vcf_path)

        self.assertEqual(alignment.sample_names, ["SAMPLE_REF", "SAMPLE_ALT", "SAMPLE_MISSING"])
        self.assertEqual(alignment.sequences["SAMPLE_REF"], "AC")
        self.assertEqual(alignment.sequences["SAMPLE_ALT"], "GT")
        self.assertEqual(alignment.sequences["SAMPLE_MISSING"], "NN")
        self.assertEqual(len(alignment.sites), 2)
        self.assertEqual(alignment.sites[0]["position"], "10")


class SnpAlignmentOutputTests(unittest.TestCase):
    def test_write_alignment_outputs_records_fasta_site_table_summary_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            input_row = SnpAlignmentInput(
                organelle="cpDNA",
                track_id="cpdna_population_sites",
                sample_count=3,
                filtered_records=2,
                filtered_vcf_path=Path("cpDNA.primary.filtered.vcf.gz"),
                filtered_vcf_index_path=Path("cpDNA.primary.filtered.vcf.gz.tbi"),
            )
            alignment = build_snp_alignment_from_literals(
                sample_names=["A", "B"],
                sequences={"A": "AC", "B": "GN"},
            )

            write_alignment_outputs(
                output_dir=output_dir,
                results=[
                    input_row.to_result(
                        alignment=alignment,
                        alignment_fasta_path=output_dir / "cpDNA.primary.snp_alignment.fa",
                        site_table_path=output_dir / "cpDNA.primary.snp_sites.tsv",
                    )
                ],
                run_label="primary",
            )

            summary = (output_dir / "primary.snp_alignment_summary.tsv").read_text()
            report = (output_dir / "primary.snp_alignment_report.md").read_text()

        self.assertIn("alignment_sites", summary)
        self.assertIn("missing_bases", summary)
        self.assertIn("# SNP Alignment", report)
        self.assertIn("SNP-only FASTA alignments", report)


def build_snp_alignment_from_literals(sample_names, sequences):
    from dudleya_organelle_alignment_pipeline.snp_alignment import SnpAlignment

    return SnpAlignment(
        sample_names=sample_names,
        sequences=sequences,
        sites=[
            {
                "site_index": "1",
                "chrom": "chloroplast",
                "position": "10",
                "ref": "A",
                "alt": "G",
            },
            {
                "site_index": "2",
                "chrom": "chloroplast",
                "position": "20",
                "ref": "T",
                "alt": "C",
            },
        ],
    )


if __name__ == "__main__":
    unittest.main()
