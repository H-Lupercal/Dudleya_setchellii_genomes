import gzip
import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.callable_consensus import (
    ConsensusInput,
    build_callable_consensus,
    write_consensus_outputs,
)


class CallableConsensusBuildTests(unittest.TestCase):
    def test_build_callable_consensus_applies_depth_variants_and_failed_site_mask(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference_path = root / "reference.fa"
            reference_path.write_text(">chloroplast\nACGTACGT\n")
            bed_path = root / "regions.bed"
            bed_path.write_text("chloroplast\t1\t7\tcallable\n")
            sample_table = root / "included_samples.tsv"
            sample_table.write_text(
                "sample_id\tdownstream_cpDNA_use\tdownstream_mtDNA_use\n"
                "S1\tinclude\tinclude\n"
                "S2\tinclude\tinclude\n"
            )
            depth_dir = root / "qc"
            depth_dir.mkdir()
            (depth_dir / "S1.depth.tsv").write_text(
                "chloroplast\t2\t3\n"
                "chloroplast\t3\t3\n"
                "chloroplast\t4\t3\n"
                "chloroplast\t6\t3\n"
                "chloroplast\t7\t3\n"
            )
            (depth_dir / "S2.depth.tsv").write_text(
                "chloroplast\t2\t3\n"
                "chloroplast\t3\t3\n"
                "chloroplast\t4\t3\n"
                "chloroplast\t5\t3\n"
                "chloroplast\t6\t3\n"
                "chloroplast\t7\t3\n"
            )
            filtered_vcf = root / "filtered.vcf.gz"
            with gzip.open(filtered_vcf, "wt") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                )
                handle.write("chloroplast\t3\t.\tG\tT\t.\tPASS\t.\tGT\t1\t0\n")
            raw_vcf = root / "raw.vcf.gz"
            with gzip.open(raw_vcf, "wt") as handle:
                handle.write("##fileformat=VCFv4.2\n")
                handle.write(
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
                )
                handle.write("chloroplast\t3\t.\tG\tT\t.\tPASS\t.\tGT\t1\t0\n")
                handle.write("chloroplast\t6\t.\tC\tA\t.\tLowQual\t.\tGT\t1\t1\n")

            alignment = build_callable_consensus(
                reference_path=reference_path,
                bed_path=bed_path,
                sample_table=sample_table,
                depth_dir=depth_dir,
                raw_vcf_path=raw_vcf,
                filtered_vcf_path=filtered_vcf,
                min_depth=1,
                organelle="cpDNA",
            )

        self.assertEqual(alignment.sample_names, ["S1", "S2"])
        self.assertEqual(alignment.consensus_length, 6)
        self.assertEqual(alignment.sequences["S1"], "CTTNNG")
        self.assertEqual(alignment.sequences["S2"], "CGTANG")
        self.assertEqual(alignment.filtered_variant_sites, 1)
        self.assertEqual(alignment.masked_failed_variant_sites, 1)
        self.assertEqual(alignment.missing_bases, 3)


class CallableConsensusOutputTests(unittest.TestCase):
    def test_write_consensus_outputs_records_summary_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            result = ConsensusInput(
                organelle="cpDNA",
                track_id="cpdna_population_sites",
                bed_path=Path("cpdna_population_sites.bed"),
                raw_records=2,
                filtered_records=1,
                raw_vcf_path=Path("raw.vcf.gz"),
                filtered_vcf_path=Path("filtered.vcf.gz"),
            ).to_result(
                sample_count=2,
                consensus_length=6,
                filtered_variant_sites=1,
                masked_failed_variant_sites=1,
                missing_bases=3,
                fasta_path=output_dir / "cpDNA.primary.callable_consensus.fa",
                site_table_path=output_dir / "cpDNA.primary.callable_sites.tsv",
            )

            write_consensus_outputs(output_dir, [result], run_label="primary", min_depth=1)

            summary = (output_dir / "primary.callable_consensus_summary.tsv").read_text()
            report = (output_dir / "primary.callable_consensus_report.md").read_text()

        self.assertIn("consensus_length", summary)
        self.assertIn("masked_failed_variant_sites", summary)
        self.assertIn("# Callable-Site Consensus Alignment", report)
        self.assertIn("full callable-site FASTA alignments", report)


if __name__ == "__main__":
    unittest.main()
