import csv
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from dudleya_organelle_alignment_pipeline.haplotype_network import (
    HaplotypeNetworkError,
    build_popart_nexus,
    build_renderer_command,
    filter_complete_case_sites,
    network_paths,
    run_haplotype_network_analysis,
    validate_sample_metadata,
    validate_renderer_outputs,
    write_network_input_fasta,
    write_network_metadata,
    write_network_site_table,
)


class CompleteCaseTests(unittest.TestCase):
    def test_filter_complete_case_sites_drops_any_non_acgt_column(self):
        records = [("S1", "ACNT"), ("S2", "ATGT"), ("S3", "ACGT")]

        filtered, kept, dropped = filter_complete_case_sites(records)

        self.assertEqual(
            filtered,
            [("S1", "ACT"), ("S2", "ATT"), ("S3", "ACT")],
        )
        self.assertEqual(kept, [0, 1, 3])
        self.assertEqual(dropped, [2])

    def test_validate_sample_metadata_rejects_mismatch(self):
        with self.assertRaisesRegex(HaplotypeNetworkError, "sample IDs"):
            validate_sample_metadata(["S1", "S2"], {"S1": {}, "S3": {}})

    def test_input_writers_record_site_filter_and_unresolved_species(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            fasta_path = output_dir / "network.fa"
            site_path = output_dir / "sites.tsv"
            metadata_path = output_dir / "metadata.tsv"
            records = [("S1", "ACT"), ("S2", "ATT")]
            source_rows = [
                {"organelle": "cpDNA", "position": str(position)}
                for position in (10, 20, 30, 40)
            ]
            metadata = {
                "S1": {
                    "species": "D. cymosa",
                    "popcode": "CY_RED",
                    "population_name": "Red",
                    "naming_profile": "main_standard",
                },
                "S2": {
                    "species": "",
                    "popcode": "",
                    "population_name": "",
                    "naming_profile": "initial_du_dash",
                },
            }

            write_network_input_fasta(fasta_path, records)
            write_network_site_table(site_path, source_rows, [0, 1, 3])
            write_network_metadata(metadata_path, records, metadata)

            fasta_text = fasta_path.read_text()
            with site_path.open(newline="") as handle:
                site_rows = list(csv.DictReader(handle, delimiter="\t"))
            with metadata_path.open(newline="") as handle:
                metadata_rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(
            [row["network_status"] for row in site_rows],
            ["retained", "retained", "dropped_missing", "retained"],
        )
        self.assertEqual(fasta_text, ">S1\nACT\n>S2\nATT\n")
        self.assertEqual(metadata_rows[1]["species_group"], "unresolved")


class ExportAndRendererContractTests(unittest.TestCase):
    def test_popart_nexus_and_renderer_command(self):
        metadata = {
            "S1": {"species": "D. cymosa"},
            "S2": {"species": ""},
        }

        nexus = build_popart_nexus(
            [("S1", "ACT"), ("S2", "ATT")],
            metadata,
        )
        command = build_renderer_command(
            Path("/bin/Rscript"),
            Path("render.R"),
            Path("in.fa"),
            Path("meta.tsv"),
            Path("out/cpDNA.primary"),
            "cpDNA",
        )

        self.assertIn("#NEXUS", nexus)
        self.assertIn("NTAX=2 NCHAR=3", nexus)
        self.assertIn("BEGIN TRAITS;", nexus)
        self.assertIn("DIMENSIONS NTRAITS=2", nexus)
        self.assertIn("TRAITLABELS D._cymosa unresolved;", nexus)
        self.assertIn("S1 1,0", nexus)
        self.assertIn("S2 0,1", nexus)
        self.assertEqual(
            command,
            [
                "/bin/Rscript",
                "render.R",
                "in.fa",
                "meta.tsv",
                "out/cpDNA.primary",
                "cpDNA",
            ],
        )

    def test_renderer_validation_rejects_missing_or_empty_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = network_paths(Path(tmp), "cpDNA", "primary")
            for path in paths.renderer_outputs:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("content\n")

            paths.edges.write_text("")
            with self.assertRaisesRegex(
                HaplotypeNetworkError,
                str(paths.edges),
            ):
                validate_renderer_outputs(paths)

            paths.edges.write_text("content\n")
            paths.svg.unlink()
            with self.assertRaisesRegex(
                HaplotypeNetworkError,
                str(paths.svg),
            ):
                validate_renderer_outputs(paths)


class PegasRendererIntegrationTests(unittest.TestCase):
    def test_renderer_writes_network_tables_and_three_figure_formats(self):
        pipeline_dir = Path(__file__).resolve().parents[1]
        repo_root = pipeline_dir.parent
        rscript = repo_root / ".tools/bioconda-env/bin/Rscript"
        renderer = pipeline_dir / "scripts/render_haplotype_network.R"

        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            fasta = output_dir / "tiny.fa"
            metadata = output_dir / "metadata.tsv"
            paths = network_paths(output_dir, "cpDNA", "tiny")
            fasta.write_text(">S1\nAAA\n>S2\nAAA\n>S3\nAAT\n>S4\nATT\n")
            metadata.write_text(
                "sample_id\tspecies_group\tspecies\tpopcode\tpopulation_name\tnaming_profile\n"
                "S1\tD. cymosa\tD. cymosa\tCY1\tOne\tmain_standard\n"
                "S2\tD. cymosa\tD. cymosa\tCY1\tOne\tmain_standard\n"
                "S3\tD. setchellii\tD. setchellii\tSE1\tTwo\tmain_standard\n"
                "S4\tD. setchellii\tD. setchellii\tSE1\tTwo\tmain_standard\n"
            )

            completed = subprocess.run(
                build_renderer_command(
                    rscript,
                    renderer,
                    fasta,
                    metadata,
                    paths.prefix,
                    "cpDNA",
                ),
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            validate_renderer_outputs(paths)
            with paths.assignments.open(newline="") as handle:
                assignments = list(csv.DictReader(handle, delimiter="\t"))
            with paths.haplotype_summary.open(newline="") as handle:
                haplotypes = list(csv.DictReader(handle, delimiter="\t"))
            with paths.edges.open(newline="") as handle:
                edges = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(len(assignments), 4)
        self.assertEqual(len(haplotypes), 3)
        self.assertGreater(len(edges), 0)


class Stage21IntegrationTests(unittest.TestCase):
    def test_stage_writes_combined_summary_commands_and_report(self):
        def renderer_stub(command, **kwargs):
            self.assertTrue(kwargs["text"])
            self.assertTrue(kwargs["capture_output"])
            self.assertFalse(kwargs["check"])
            prefix = Path(command[4])
            organelle = command[5]
            Path(f"{prefix}.haplotype_assignments.tsv").write_text(
                "sample_id\torganelle\thaplotype_id\tspecies_group\tpopcode\n"
                f"S1\t{organelle}\tH001\tD. cymosa\tCY1\n"
                f"S2\t{organelle}\tH002\tD. setchellii\tSE1\n"
                f"S3\t{organelle}\tH001\tunresolved\t\n"
            )
            Path(f"{prefix}.haplotype_summary.tsv").write_text(
                "organelle\thaplotype_id\tsample_count\tspecies_group_count\tpopcode_count\n"
                f"{organelle}\tH001\t2\t2\t1\n"
                f"{organelle}\tH002\t1\t1\t1\n"
            )
            Path(f"{prefix}.haplotype_network_edges.tsv").write_text(
                "organelle\tfrom_haplotype\tto_haplotype\tmutation_steps\talternative_link\n"
                f"{organelle}\tH002\tH001\t1\tFALSE\n"
            )
            Path(f"{prefix}.haplotype_network_layout.tsv").write_text(
                "organelle\thaplotype_id\tx\ty\n"
                f"{organelle}\tH001\t-1\t0\n"
                f"{organelle}\tH002\t1\t0\n"
            )
            Path(f"{prefix}.haplotype_network_renderer_summary.tsv").write_text(
                "organelle\tsample_count\thaplotype_count\tedge_count\tspecies_group_count\n"
                f"{organelle}\t3\t2\t1\t3\n"
            )
            for suffix in ("png", "pdf", "svg"):
                Path(f"{prefix}.haplotype_network.{suffix}").write_bytes(b"figure")
            return SimpleNamespace(returncode=0, stdout="", stderr="")

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            snp_dir = root / "stage10"
            output_dir = root / "stage21"
            snp_dir.mkdir()
            metadata_path = root / "included_samples.tsv"
            metadata_path.write_text(
                "sample_id\tspecies\tpopcode\tpopulation_name\tnaming_profile\n"
                "S1\tD. cymosa\tCY1\tOne\tmain_standard\n"
                "S2\tD. setchellii\tSE1\tTwo\tmain_standard\n"
                "S3\t\t\t\tinitial_du_dash\n"
            )
            inputs = {
                "cpDNA": [("S1", "ACNT"), ("S2", "ATGT"), ("S3", "ACGT")],
                "mtDNA": [("S1", "AAAN"), ("S2", "AATT"), ("S3", "AACT")],
            }
            summary_rows = []
            for organelle, records in inputs.items():
                fasta_path = snp_dir / f"{organelle}.primary.snp_alignment.fa"
                site_path = snp_dir / f"{organelle}.primary.snp_sites.tsv"
                fasta_path.write_text(
                    "".join(
                        f">{sample_id}\n{sequence}\n"
                        for sample_id, sequence in records
                    )
                )
                site_path.write_text(
                    "site_index\tchrom\tposition\tref\talt\n"
                    + "".join(
                        f"{index}\t{organelle}\t{index * 10}\tA\tT\n"
                        for index in range(1, 5)
                    )
                )
                summary_rows.append(
                    f"{organelle}\t{organelle.lower()}_track\t3\t4\t4\t1\tunused.vcf.gz\t{fasta_path}\t{site_path}\n"
                )
            (snp_dir / "primary.snp_alignment_summary.tsv").write_text(
                "organelle\ttrack_id\tsample_count\tfiltered_records\talignment_sites\tmissing_bases\tfiltered_vcf_path\talignment_fasta_path\tsite_table_path\n"
                + "".join(summary_rows)
            )

            results = run_haplotype_network_analysis(
                snp_alignment_dir=snp_dir,
                metadata_path=metadata_path,
                output_dir=output_dir,
                run_label="primary",
                rscript=Path("Rscript"),
                renderer_path=Path("renderer.R"),
                runner=renderer_stub,
            )
            summary = (output_dir / "primary.haplotype_network_summary.tsv").read_text()
            commands = (output_dir / "primary.haplotype_network_commands.tsv").read_text()
            report = (output_dir / "primary.haplotype_network_report.md").read_text()

        self.assertEqual([result.organelle for result in results], ["cpDNA", "mtDNA"])
        self.assertIn("retained_site_count", summary)
        self.assertIn("Rscript renderer.R", commands)
        self.assertIn("cpDNA", report)
        self.assertIn("mtDNA", report)
        self.assertIn("pegas::haploNet", report)
        self.assertIn("complete-case filtering", report)
        self.assertIn("not ancestry proportions", report)
        self.assertIn("primary links only", report)
        self.assertIn("alternative links", report)


if __name__ == "__main__":
    unittest.main()
