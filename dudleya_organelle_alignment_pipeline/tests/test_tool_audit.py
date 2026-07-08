import tempfile
import unittest
from pathlib import Path

from dudleya_organelle_alignment_pipeline.tool_audit import (
    ToolSpec,
    check_tool,
    summarize_audit,
    write_tool_audit_outputs,
)


class ToolAuditCheckTests(unittest.TestCase):
    def test_check_tool_records_found_tool_path_and_version(self):
        spec = ToolSpec(
            tool_id="samtools",
            executables=("samtools",),
            necessity="required_current",
            required_for="BAM processing",
            version_args=("--version",),
        )

        result = check_tool(
            spec,
            resolver=lambda executable: "/env/bin/samtools" if executable == "samtools" else None,
            runner=lambda command: "samtools 1.23.1\nUsing htslib",
        )

        self.assertEqual(result.status, "FOUND")
        self.assertEqual(result.executable, "samtools")
        self.assertEqual(result.path, "/env/bin/samtools")
        self.assertEqual(result.version, "samtools 1.23.1")

    def test_check_tool_records_missing_required_tool(self):
        spec = ToolSpec(
            tool_id="admixture",
            executables=("admixture",),
            necessity="required_remaining",
            required_for="structure/admixture-style clustering",
            version_args=("--version",),
        )

        result = check_tool(spec, resolver=lambda executable: None, runner=lambda command: "")

        self.assertEqual(result.status, "MISSING")
        self.assertEqual(result.path, "")
        self.assertIn("Install before", result.note)

    def test_check_tool_marks_found_executable_missing_when_package_check_fails(self):
        spec = ToolSpec(
            tool_id="python_matplotlib",
            executables=("python3",),
            necessity="required_remaining",
            required_for="PCA and tree figure rendering",
            version_args=("-c", "import matplotlib; print(matplotlib.__version__)"),
        )

        def failing_runner(command):
            raise RuntimeError("No module named matplotlib")

        result = check_tool(
            spec,
            resolver=lambda executable: "/env/bin/python3",
            runner=failing_runner,
        )

        self.assertEqual(result.status, "MISSING")
        self.assertEqual(result.path, "/env/bin/python3")
        self.assertIn("version/import check failed", result.note)

    def test_visualization_dependencies_are_in_default_audit_specs(self):
        from dudleya_organelle_alignment_pipeline.tool_audit import TOOL_SPECS

        tool_ids = {spec.tool_id for spec in TOOL_SPECS}

        self.assertIn("python_matplotlib", tool_ids)
        self.assertIn("python_pandas", tool_ids)
        self.assertIn("python_sklearn", tool_ids)
        self.assertIn("python_biopython", tool_ids)
        self.assertIn("r_ggplot2", tool_ids)
        self.assertIn("r_ape", tool_ids)

    def test_summarize_audit_flags_missing_required_remaining_tools(self):
        specs = [
            ToolSpec(
                tool_id="iqtree",
                executables=("iqtree",),
                necessity="required_current",
                required_for="ML trees",
                version_args=("--version",),
            ),
            ToolSpec(
                tool_id="admixture",
                executables=("admixture",),
                necessity="required_remaining",
                required_for="structure/admixture-style clustering",
                version_args=("--version",),
            ),
        ]
        results = [
            check_tool(specs[0], resolver=lambda executable: "/bin/iqtree", runner=lambda command: "IQ-TREE 3"),
            check_tool(specs[1], resolver=lambda executable: None, runner=lambda command: ""),
        ]

        summary = summarize_audit(results)

        self.assertFalse(summary.ready_for_remaining_goal)
        self.assertEqual(summary.missing_required_remaining, ["admixture"])


class ToolAuditOutputTests(unittest.TestCase):
    def test_write_tool_audit_outputs_records_tsv_and_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            spec = ToolSpec(
                tool_id="bcftools",
                executables=("bcftools",),
                necessity="required_current",
                required_for="haploid variant calling",
                version_args=("--version",),
            )
            result = check_tool(
                spec,
                resolver=lambda executable: "/env/bin/bcftools",
                runner=lambda command: "bcftools 1.23.1",
            )

            write_tool_audit_outputs(output_dir, [result], audit_label="primary")

            tsv = (output_dir / "primary.tool_audit.tsv").read_text()
            report = (output_dir / "primary.tool_audit_report.md").read_text()

        self.assertIn("tool_id", tsv)
        self.assertIn("bcftools", tsv)
        self.assertIn("# Bioinformatics Tool Audit", report)
        self.assertIn("bcftools", report)


if __name__ == "__main__":
    unittest.main()
