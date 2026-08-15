import importlib.util
import unittest
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "run_phylogenetic_tree_10000.py"
)


class TenThousandBootstrapRunnerTests(unittest.TestCase):
    def test_runner_uses_10000_replicates_14_threads_and_separate_output(self):
        spec = importlib.util.spec_from_file_location(
            "run_phylogenetic_tree_10000",
            SCRIPT_PATH,
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        arguments = module.build_run_arguments()

        self.assertEqual(
            arguments[arguments.index("--bootstrap-replicates") + 1],
            "10000",
        )
        self.assertEqual(arguments[arguments.index("--threads") + 1], "14")
        self.assertEqual(
            arguments[arguments.index("--output-dir") + 1],
            "dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree_10000",
        )


if __name__ == "__main__":
    unittest.main()
