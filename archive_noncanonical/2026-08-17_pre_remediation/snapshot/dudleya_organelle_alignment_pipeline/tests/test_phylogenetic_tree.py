import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.phylogenetic_tree import (
    TreeInput,
    build_iqtree_command,
    tree_output_prefix,
    write_tree_outputs,
)


class PhylogeneticTreeCommandTests(unittest.TestCase):
    def test_build_iqtree_command_uses_callable_alignment_and_fast_ml_options(self):
        command = build_iqtree_command(
            iqtree_executable="iqtree",
            alignment_path=Path("cpDNA.primary.callable_consensus.fa"),
            prefix=Path("results/cpDNA.primary.iqtree_ml"),
            model="GTR+F+G4",
            threads=4,
            fast=True,
            bootstrap_replicates=0,
        )

        self.assertEqual(command[0], "iqtree")
        self.assertIn("-s", command)
        self.assertIn("cpDNA.primary.callable_consensus.fa", command)
        self.assertIn("--seqtype", command)
        self.assertIn("DNA", command)
        self.assertIn("-m", command)
        self.assertIn("GTR+F+G4", command)
        self.assertIn("--fast", command)
        self.assertIn("--redo", command)

    def test_build_iqtree_command_can_request_ultrafast_bootstrap_support(self):
        command = build_iqtree_command(
            iqtree_executable="iqtree",
            alignment_path=Path("mtDNA.primary.callable_consensus.fa"),
            prefix=Path("results/mtDNA.primary.iqtree_ml_bootstrap"),
            model="GTR+F+G4",
            threads=4,
            fast=False,
            bootstrap_replicates=1000,
        )

        self.assertIn("-B", command)
        self.assertIn("1000", command)
        self.assertIn("--bnni", command)
        self.assertNotIn("--fast", command)

    def test_tree_output_prefix_includes_organelle_run_label_and_method(self):
        tree_input = TreeInput(
            organelle="mtDNA",
            track_id="mtdna_high_confidence_unique",
            sample_count=275,
            alignment_sites=44930,
            missing_bases=31313,
            alignment_fasta_path=Path("mtDNA.primary.callable_consensus.fa"),
        )

        prefix = tree_output_prefix(tree_input, Path("trees"), run_label="primary")

        self.assertEqual(prefix, Path("trees/mtDNA.primary.iqtree_ml"))


class PhylogeneticTreeOutputTests(unittest.TestCase):
    def test_write_tree_outputs_records_summary_commands_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            result = TreeInput(
                organelle="cpDNA",
                track_id="cpdna_population_sites",
                sample_count=2,
                alignment_sites=8,
                missing_bases=0,
                alignment_fasta_path=Path("cpDNA.primary.callable_consensus.fa"),
            ).to_result(
                model="GTR+F+G4",
                method="iqtree_ml_fast",
                tree_prefix=output_dir / "cpDNA.primary.iqtree_ml",
                treefile_path=output_dir / "cpDNA.primary.iqtree_ml.treefile",
                log_path=output_dir / "cpDNA.primary.iqtree_ml.log",
                iqtree_report_path=output_dir / "cpDNA.primary.iqtree_ml.iqtree",
            )

            write_tree_outputs(
                output_dir=output_dir,
                results=[result],
                command_rows=[
                    {
                        "organelle": "cpDNA",
                        "method": "iqtree_ml_fast",
                        "command": "iqtree -s cpDNA.primary.callable_consensus.fa",
                    }
                ],
                run_label="primary",
            )

            summary = (output_dir / "primary.phylogenetic_tree_summary.tsv").read_text()
            report = (output_dir / "primary.phylogenetic_tree_report.md").read_text()
            commands = (output_dir / "primary.phylogenetic_tree_commands.tsv").read_text()

        self.assertIn("treefile_path", summary)
        self.assertIn("# Phylogenetic Trees", report)
        self.assertIn("IQ-TREE maximum-likelihood", report)
        self.assertIn("iqtree -s", commands)


if __name__ == "__main__":
    unittest.main()
