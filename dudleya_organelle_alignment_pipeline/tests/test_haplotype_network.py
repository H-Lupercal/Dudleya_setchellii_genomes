import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.haplotype_network import (
    HaplotypeNetworkError,
    build_popart_nexus,
    build_renderer_command,
    filter_complete_case_sites,
    network_paths,
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
        self.assertIn("S1 D._cymosa", nexus)
        self.assertIn("S2 unresolved", nexus)
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


if __name__ == "__main__":
    unittest.main()
