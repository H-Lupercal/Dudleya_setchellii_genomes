import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.tree_visualization import (
    TreeFigureInput,
    compute_tree_figure_size,
    read_tree_figure_inputs,
    render_tree_figure,
    write_tree_visualization_outputs,
)


class TreeVisualizationInputTests(unittest.TestCase):
    def test_read_tree_figure_inputs_uses_phylogenetic_tree_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            tree_dir = Path(tmp)
            tree_path = tree_dir / "cpDNA.primary.iqtree_ml.treefile"
            tree_path.write_text("(DU-1:0.1,DU-2:0.2);\n")
            (tree_dir / "primary.phylogenetic_tree_summary.tsv").write_text(
                "organelle\ttrack_id\tsample_count\talignment_sites\tmissing_bases\t"
                "method\tmodel\talignment_fasta_path\ttreefile_path\tlog_path\t"
                "iqtree_report_path\n"
                f"cpDNA\tcpdna_population_sites\t2\t8\t0\tiqtree_ml_fast\t"
                f"GTR+F+G4\tcp.fa\t{tree_path}\tcp.log\tcp.iqtree\n"
            )

            inputs = read_tree_figure_inputs(tree_dir, run_label="primary")

        self.assertEqual(len(inputs), 1)
        self.assertEqual(inputs[0].organelle, "cpDNA")
        self.assertEqual(inputs[0].sample_count, 2)
        self.assertEqual(inputs[0].treefile_path, tree_path)


class TreeVisualizationRenderingTests(unittest.TestCase):
    def test_compute_tree_figure_size_grows_for_many_samples(self):
        small = compute_tree_figure_size(sample_count=5)
        large = compute_tree_figure_size(sample_count=275)

        self.assertGreater(large[1], small[1])
        self.assertGreaterEqual(large[0], 10)

    def test_render_tree_figure_writes_png_pdf_and_svg(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            tree_path = output_dir / "mtDNA.primary.iqtree_ml.treefile"
            tree_path.write_text("(DU-1:0.1,(DU-2:0.2,DU-3:0.3):0.4);\n")
            figure_input = TreeFigureInput(
                organelle="mtDNA",
                track_id="mtdna_high_confidence_unique",
                sample_count=3,
                alignment_sites=12,
                model="GTR+F+G4",
                method="iqtree_ml_fast",
                treefile_path=tree_path,
            )

            result = render_tree_figure(figure_input, output_dir, run_label="primary")

            self.assertTrue(result.png_path.exists())
            self.assertTrue(result.pdf_path.exists())
            self.assertTrue(result.svg_path.exists())
            self.assertGreater(result.tip_count, 0)

    def test_write_tree_visualization_outputs_records_summary_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            tree_path = output_dir / "cpDNA.primary.iqtree_ml.treefile"
            tree_path.write_text("(DU-1:0.1,DU-2:0.2);\n")
            figure_input = TreeFigureInput(
                organelle="cpDNA",
                track_id="cpdna_population_sites",
                sample_count=2,
                alignment_sites=8,
                model="GTR+F+G4",
                method="iqtree_ml_fast",
                treefile_path=tree_path,
            )
            result = render_tree_figure(figure_input, output_dir, run_label="primary")

            write_tree_visualization_outputs(output_dir, [result], run_label="primary")

            summary = (output_dir / "primary.tree_visualization_summary.tsv").read_text()
            report = (output_dir / "primary.tree_visualization_report.md").read_text()

        self.assertIn("png_path", summary)
        self.assertIn("# Step 14 Tree Visualizations", report)
        self.assertIn("cpDNA", report)


if __name__ == "__main__":
    unittest.main()
