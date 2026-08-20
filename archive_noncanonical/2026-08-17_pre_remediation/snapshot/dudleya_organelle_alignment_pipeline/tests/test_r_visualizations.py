import csv
import subprocess
import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.r_visualizations import (
    FigureJob,
    build_arg_parser,
    build_renderer_command,
    discover_figure_jobs,
    figure_outputs,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
PIPELINE_DIR = REPO_ROOT / "dudleya_organelle_alignment_pipeline"
RSCRIPT = REPO_ROOT / ".tools/bioconda-env/bin/Rscript"


class RVisualizationOrchestrationTests(unittest.TestCase):
    def test_figure_outputs_are_additive_and_include_three_formats(self):
        outputs = figure_outputs(Path("results/cpDNA.primary.pca"), "r_ggplot")

        self.assertEqual(
            outputs,
            (
                Path("results/cpDNA.primary.pca.r_ggplot.png"),
                Path("results/cpDNA.primary.pca.r_ggplot.pdf"),
                Path("results/cpDNA.primary.pca.r_ggplot.svg"),
            ),
        )

    def test_build_renderer_command_uses_source_data_and_additive_prefix(self):
        job = FigureJob(
            family="pca",
            organelle="cpDNA",
            stage="15_pca",
            renderer_path=Path("scripts/render_pca_ggplot.R"),
            renderer_suffix="r_ggplot",
            source_paths=(Path("coordinates.tsv"), Path("variance.tsv")),
            output_prefix=Path("results/cpDNA.primary.pca"),
            extra_args=("cpDNA",),
        )

        command = build_renderer_command(Path("Rscript"), job)

        self.assertEqual(
            command,
            [
                "Rscript",
                "scripts/render_pca_ggplot.R",
                "coordinates.tsv",
                "variance.tsv",
                "results/cpDNA.primary.pca.r_ggplot",
                "cpDNA",
            ],
        )
        self.assertNotIn("results/cpDNA.primary.pca.png", command)

    def test_discover_jobs_covers_every_existing_non_r_figure(self):
        jobs = discover_figure_jobs(
            pipeline_dir=PIPELINE_DIR,
            run_label="primary",
        )

        families = [job.family for job in jobs]
        self.assertEqual(families.count("pca"), 2)
        self.assertEqual(families.count("admixture_structure"), 4)
        self.assertEqual(families.count("admixture_cv"), 4)
        self.assertEqual(families.count("tree"), 4)
        self.assertEqual(len(jobs), 14)
        self.assertTrue(all(job.renderer_suffix in {"r_ggplot", "r_ggtree"} for job in jobs))

    def test_cli_accepts_rscript_run_label_and_stage_selection(self):
        args = build_arg_parser().parse_args(
            [
                "--rscript",
                "/opt/Rscript",
                "--run-label",
                "review",
                "--stages",
                "15_pca",
                "20_bootstrap_tree_visualization",
            ]
        )

        self.assertEqual(args.rscript, Path("/opt/Rscript"))
        self.assertEqual(args.run_label, "review")
        self.assertEqual(
            args.stages,
            ["15_pca", "20_bootstrap_tree_visualization"],
        )

    def test_runner_script_can_be_executed_directly_from_repo_root(self):
        completed = subprocess.run(
            [
                "python3",
                str(PIPELINE_DIR / "scripts/run_r_visualizations.py"),
                "--help",
            ],
            cwd=REPO_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )

        self.assertEqual(completed.returncode, 0, completed.stdout)
        self.assertIn("additive R alternatives", completed.stdout)


@unittest.skipUnless(RSCRIPT.exists(), "pipeline R environment is unavailable")
class RVisualizationRendererSmokeTests(unittest.TestCase):
    def run_renderer(self, script_name: str, args: list[Path | str]) -> None:
        command = [str(RSCRIPT), str(PIPELINE_DIR / "scripts" / script_name)]
        command.extend(str(arg) for arg in args)
        completed = subprocess.run(
            command,
            cwd=REPO_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout)

    def assert_figure_set(self, prefix: Path) -> None:
        for suffix in ("png", "pdf", "svg"):
            path = Path(f"{prefix}.{suffix}")
            self.assertTrue(path.exists(), path)
            self.assertGreater(path.stat().st_size, 100, path)

    def test_pca_renderer_writes_all_formats_with_species_legend(self):
        with tempfile.TemporaryDirectory() as temp_name:
            temp_dir = Path(temp_name)
            coordinates = temp_dir / "coordinates.tsv"
            variance = temp_dir / "variance.tsv"
            prefix = temp_dir / "pca.r_ggplot"
            coordinates.write_text(
                "sample_id\torganelle\tpc1\tpc2\tspecies\tpopcode\tpopulation_name\t"
                "naming_profile\tplot_group\n"
                "s1\tcpDNA\t-1\t0.5\tD. setchellii\tA\tAlpha\tx\tg1\n"
                "s2\tcpDNA\t1\t-0.5\tD. cymosa\tB\tBeta\tx\tg2\n"
                "s3\tcpDNA\t0\t0\t\tC\tGamma\tx\tg3\n"
            )
            variance.write_text(
                "organelle\tcomponent\texplained_variance_ratio\tretained_sites\n"
                "cpDNA\tPC1\t0.42\t20\n"
                "cpDNA\tPC2\t0.21\t20\n"
            )

            self.run_renderer(
                "render_pca_ggplot.R",
                [coordinates, variance, prefix, "cpDNA"],
            )

            self.assert_figure_set(prefix)

    def test_admixture_renderer_writes_structure_and_cv_formats(self):
        with tempfile.TemporaryDirectory() as temp_name:
            temp_dir = Path(temp_name)
            q_table = temp_dir / "q.tsv"
            summary = temp_dir / "summary.tsv"
            structure_prefix = temp_dir / "structure.r_ggplot"
            cv_prefix = temp_dir / "cv.r_ggplot"
            q_table.write_text(
                "sample_id\torganelle\tbest_k\tspecies\tpopcode\tpopulation_name\t"
                "plot_group\tcluster_1\tcluster_2\n"
                "s1\tcpDNA\t2\tD. setchellii\tA\tAlpha\tg1\t0.9\t0.1\n"
                "s2\tcpDNA\t2\tD. setchellii\tA\tAlpha\tg1\t0.8\t0.2\n"
                "s3\tcpDNA\t2\tD. cymosa\tB\tBeta\tg2\t0.2\t0.8\n"
            )
            with summary.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["organelle", "k", "replicate", "cv_error", "is_best_mean_k"],
                    delimiter="\t",
                )
                writer.writeheader()
                for k, values in ((1, (0.30, 0.31)), (2, (0.20, 0.22)), (3, (0.24, 0.25))):
                    for replicate, value in enumerate(values, start=1):
                        writer.writerow(
                            {
                                "organelle": "cpDNA",
                                "k": k,
                                "replicate": replicate,
                                "cv_error": value,
                                "is_best_mean_k": "yes" if k == 2 else "no",
                            }
                        )

            self.run_renderer(
                "render_admixture_ggplot.R",
                ["structure", q_table, structure_prefix, "cpDNA"],
            )
            self.run_renderer(
                "render_admixture_ggplot.R",
                ["cv", summary, cv_prefix, "cpDNA"],
            )

            self.assert_figure_set(structure_prefix)
            self.assert_figure_set(cv_prefix)

    def test_tree_renderer_writes_initial_and_bootstrap_formats(self):
        with tempfile.TemporaryDirectory() as temp_name:
            temp_dir = Path(temp_name)
            tree = temp_dir / "tree.nwk"
            metadata = temp_dir / "metadata.tsv"
            initial_prefix = temp_dir / "initial.r_ggtree"
            bootstrap_prefix = temp_dir / "bootstrap.r_ggtree"
            tree.write_text("((s1:0.1,s2:0.1)95:0.2,s3:0.3);\n")
            metadata.write_text(
                "sample_id\tspecies\tpopcode\n"
                "s1\tD. setchellii\tA\n"
                "s2\tD. cymosa\tB\n"
                "s3\t\tC\n"
            )

            self.run_renderer(
                "render_tree_ggtree.R",
                [tree, metadata, initial_prefix, "cpDNA", "initial", "0"],
            )
            self.run_renderer(
                "render_tree_ggtree.R",
                [tree, metadata, bootstrap_prefix, "cpDNA", "bootstrap", "1000"],
            )

            self.assert_figure_set(initial_prefix)
            self.assert_figure_set(bootstrap_prefix)


if __name__ == "__main__":
    unittest.main()
