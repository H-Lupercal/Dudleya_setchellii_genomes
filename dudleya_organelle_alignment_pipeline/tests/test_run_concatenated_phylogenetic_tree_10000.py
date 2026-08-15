import importlib.util
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "run_concatenated_phylogenetic_tree_10000.py"
)


class ConcatenatedTenThousandBootstrapRunnerTests(unittest.TestCase):
    def test_runner_uses_combined_input_10000_replicates_and_14_threads(self):
        spec = importlib.util.spec_from_file_location(
            "run_concatenated_phylogenetic_tree_10000",
            SCRIPT_PATH,
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        arguments = module.build_run_arguments()

        def value(option: str) -> str:
            return arguments[arguments.index(option) + 1]

        self.assertEqual(
            value("--consensus-dir"),
            "dudleya_organelle_alignment_pipeline/results/22_concatenated_consensus",
        )
        self.assertEqual(
            value("--output-dir"),
            "dudleya_organelle_alignment_pipeline/results/23_concatenated_bootstrap_phylogenetic_tree_10000",
        )
        self.assertEqual(value("--bootstrap-replicates"), "10000")
        self.assertEqual(value("--threads"), "14")
        self.assertEqual(value("--run-label"), "primary")


if __name__ == "__main__":
    unittest.main()
