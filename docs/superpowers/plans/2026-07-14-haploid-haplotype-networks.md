# Haploid Organelle Haplotype Networks Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a reproducible `pegas`-based Stage 21 that generates cpDNA and mtDNA haplotype-network tables, PopART inputs, reports, and PNG/PDF/SVG figures without removing any existing analysis.

**Architecture:** A Python module validates Stage 10/07 inputs, removes every SNP column containing a non-ACGT state, writes intermediate files, invokes one focused R renderer, validates its products, and writes the combined report. The R renderer uses `ape` and `pegas::haploNet` for domain-specific haplotype inference and plotting. Tests cover Python preparation/orchestration plus a small R integration run.

**Tech Stack:** Python 3.11+, `unittest`, Rscript, R `ape`, R `pegas`, TSV/FASTA/NEXUS, base R graphics.

---

## File Map

- Create `dudleya_organelle_alignment_pipeline/haplotype_network.py`: validation, preprocessing, R orchestration, summary, report, CLI.
- Create `dudleya_organelle_alignment_pipeline/scripts/render_haplotype_network.R`: `pegas` network construction, tables, figures.
- Create `dudleya_organelle_alignment_pipeline/scripts/run_haplotype_network.py`: runner wrapper.
- Create `dudleya_organelle_alignment_pipeline/tests/test_haplotype_network.py`: unit and integration tests.
- Modify `tool_audit.py`, `tests/test_tool_audit.py`, and `environment.yml`: declare and audit `r_pegas`.
- Modify `PROCESS.md`, pipeline `README.md`, integrated report, and final-deliverables manifest.
- Create `results/21_haplotype_network/` and regenerate the Stage 13 tool audit.

## Task 1: Declare and audit `pegas`

**Files:**
- Modify: `dudleya_organelle_alignment_pipeline/tests/test_tool_audit.py`
- Modify: `dudleya_organelle_alignment_pipeline/tool_audit.py`
- Modify: `dudleya_organelle_alignment_pipeline/environment.yml`

- [ ] **Step 1: Write the failing test**

Add to `test_visualization_dependencies_are_in_default_audit_specs`:

```python
self.assertIn("r_pegas", tool_ids)
```

- [ ] **Step 2: Verify RED**

Run `python3 -m unittest dudleya_organelle_alignment_pipeline.tests.test_tool_audit.ToolAuditCheckTests.test_visualization_dependencies_are_in_default_audit_specs -v`.

Expected: failure because `r_pegas` is absent.

- [ ] **Step 3: Implement the audit entry and environment dependency**

Add after `r_ape`:

```python
ToolSpec(
    "r_pegas",
    ("Rscript",),
    "required_current",
    "haploid cpDNA/mtDNA haplotype networks",
    ("-e", "cat(as.character(packageVersion('pegas')))"),
),
```

Add to `environment.yml`:

```yaml
  - r-pegas>=1.3
```

- [ ] **Step 4: Verify GREEN and commit**

Run the targeted test, then commit the three files with message `build: declare pegas haplotype network dependency`.

## Task 2: Build complete-case network inputs

**Files:**
- Create: `dudleya_organelle_alignment_pipeline/tests/test_haplotype_network.py`
- Create: `dudleya_organelle_alignment_pipeline/haplotype_network.py`

- [ ] **Step 1: Write failing filtering and sample-identity tests**

```python
import unittest

from dudleya_organelle_alignment_pipeline.haplotype_network import (
    HaplotypeNetworkError,
    filter_complete_case_sites,
    validate_sample_metadata,
)


class CompleteCaseTests(unittest.TestCase):
    def test_filter_complete_case_sites_drops_any_non_acgt_column(self):
        records = [("S1", "ACNT"), ("S2", "ATGT"), ("S3", "ACGT")]
        filtered, kept, dropped = filter_complete_case_sites(records)
        self.assertEqual(filtered, [("S1", "ACT"), ("S2", "ATT"), ("S3", "ACT")])
        self.assertEqual(kept, [0, 1, 3])
        self.assertEqual(dropped, [2])

    def test_validate_sample_metadata_rejects_mismatch(self):
        with self.assertRaisesRegex(HaplotypeNetworkError, "sample IDs"):
            validate_sample_metadata(["S1", "S2"], {"S1": {}, "S3": {}})
```

- [ ] **Step 2: Verify RED**

Run the new test module. Expected: import failure because the production module does not exist.

- [ ] **Step 3: Implement the minimal preparation API**

```python
from __future__ import annotations

from pathlib import Path

BASES = frozenset("ACGT")
DEFAULT_SNP_ALIGNMENT_DIR = Path("dudleya_organelle_alignment_pipeline/results/10_snp_alignment")
DEFAULT_METADATA_PATH = Path("dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/included_samples.tsv")
DEFAULT_OUTPUT_DIR = Path("dudleya_organelle_alignment_pipeline/results/21_haplotype_network")
DEFAULT_RENDERER_PATH = Path("dudleya_organelle_alignment_pipeline/scripts/render_haplotype_network.R")


class HaplotypeNetworkError(RuntimeError):
    pass


def filter_complete_case_sites(records):
    if not records:
        raise HaplotypeNetworkError("No FASTA records supplied")
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise HaplotypeNetworkError("FASTA records have inconsistent lengths")
    site_count = lengths.pop()
    kept = [i for i in range(site_count) if all(sequence[i] in BASES for _, sequence in records)]
    keep_set = set(kept)
    dropped = [i for i in range(site_count) if i not in keep_set]
    if not kept:
        raise HaplotypeNetworkError("No complete-case SNP sites remain")
    filtered = [(sample_id, "".join(sequence[i] for i in kept)) for sample_id, sequence in records]
    if len({sequence for _, sequence in filtered}) < 2:
        raise HaplotypeNetworkError("Fewer than two haplotypes remain")
    return filtered, kept, dropped


def validate_sample_metadata(sample_ids, metadata):
    if sample_ids != list(metadata):
        raise HaplotypeNetworkError("FASTA and metadata sample IDs or order differ")
```

- [ ] **Step 4: Verify GREEN**

Run the new test module. Expected: 2 tests pass.

- [ ] **Step 5: Add failing writer tests**

Create four source site rows, call the planned writers, and assert that positions 10, 20, and 40 are retained, position 30 is `dropped_missing`, the FASTA contains every sample in original order, and the metadata table has `species_group=unresolved` when species is blank.

- [ ] **Step 6: Implement input writers**

Implement `read_tsv`, `write_tsv`, `read_fasta`, `write_network_input_fasta`, `write_network_site_table`, and `write_network_metadata`. Follow existing module patterns. Site output fields are the original Stage 10 fields plus `source_alignment_index_0based` and `network_status`.

- [ ] **Step 7: Run tests and commit**

Run the Stage 21 tests and commit module plus tests with message `feat: prepare complete-case haplotype network inputs`.

## Task 3: Add PopART export and R orchestration

**Files:**
- Modify: `dudleya_organelle_alignment_pipeline/tests/test_haplotype_network.py`
- Modify: `dudleya_organelle_alignment_pipeline/haplotype_network.py`

- [ ] **Step 1: Write failing NEXUS and command tests**

```python
nexus = build_popart_nexus([("S1", "ACT"), ("S2", "ATT")], metadata)
self.assertIn("#NEXUS", nexus)
self.assertIn("NTAX=2 NCHAR=3", nexus)
self.assertIn("BEGIN TRAITS;", nexus)

command = build_renderer_command(Path("/bin/Rscript"), Path("render.R"), Path("in.fa"), Path("meta.tsv"), Path("out/cpDNA.primary"), "cpDNA")
self.assertEqual(command, ["/bin/Rscript", "render.R", "in.fa", "meta.tsv", "out/cpDNA.primary", "cpDNA"])
```

- [ ] **Step 2: Verify RED, then implement**

Implement a nucleotide `DATA` block and categorical species `TRAITS` block. Implement:

```python
def build_renderer_command(rscript, renderer, fasta, metadata, prefix, organelle):
    return [str(rscript), str(renderer), str(fasta), str(metadata), str(prefix), organelle]
```

- [ ] **Step 3: Add failing renderer-output validation tests**

Define expected assignments, haplotype summary, edge, layout, renderer-summary, PNG, PDF, and SVG paths. Assert a missing or empty file raises `HaplotypeNetworkError` naming that path.

- [ ] **Step 4: Implement `NetworkPaths` and validation, verify GREEN, commit**

Every required renderer file must exist and have nonzero size. Commit with message `feat: orchestrate haplotype network rendering`.

## Task 4: Implement and test the `pegas` renderer

**Files:**
- Create: `dudleya_organelle_alignment_pipeline/scripts/render_haplotype_network.R`
- Modify: `dudleya_organelle_alignment_pipeline/tests/test_haplotype_network.py`

- [ ] **Step 1: Install the declared package**

Run:

```bash
conda env update --prefix .tools/bioconda-env --file dudleya_organelle_alignment_pipeline/environment.yml
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" Rscript -e "cat(as.character(packageVersion('pegas')))"
```

Expected: `pegas` 1.3 or newer.

- [ ] **Step 2: Write and run a failing R integration test**

Create temporary sequences `AAA, AAA, AAT, ATT`, two species groups, invoke the renderer, and assert four assignments, three haplotypes, nonempty edges, and nonempty PNG/PDF/SVG. Do not skip when `pegas` is absent. Expected RED: renderer script missing.

- [ ] **Step 3: Implement the R renderer**

The entry point is:

```r
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: render_haplotype_network.R FASTA METADATA PREFIX ORGANELLE")
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(pegas))
fasta_path <- args[[1]]
metadata_path <- args[[2]]
prefix <- args[[3]]
organelle <- args[[4]]
dna <- read.dna(fasta_path, format = "fasta")
haps <- haplotype(dna)
net <- haploNet(haps)
```

Map `attr(haps, "index")` back to sample names. Build a fixed six-color species-frequency matrix. Set `set.seed(20260714)`, use `show.mutation=3`, scale node area by sample count, suppress an automatic legend, record returned node coordinates, and reuse those coordinates for PNG/PDF/SVG.

Write these columns:

```text
assignments: sample_id, organelle, haplotype_id, species_group, popcode
haplotype_summary: organelle, haplotype_id, sample_count, species_group_count, popcode_count
edges: organelle, from_haplotype, to_haplotype, mutation_steps, alternative_link
layout: organelle, haplotype_id, x, y
renderer_summary: organelle, sample_count, haplotype_count, edge_count, species_group_count
```

- [ ] **Step 4: Verify GREEN and commit**

Run the integration test with the pipeline environment. Inspect the tiny figures. Commit with message `feat: render pegas haplotype networks`.

## Task 5: Complete Stage 21 and its runner

**Files:**
- Modify: `dudleya_organelle_alignment_pipeline/haplotype_network.py`
- Modify: `dudleya_organelle_alignment_pipeline/tests/test_haplotype_network.py`
- Create: `dudleya_organelle_alignment_pipeline/scripts/run_haplotype_network.py`

- [ ] **Step 1: Write a failing end-to-end stage test**

Use miniature cpDNA/mtDNA Stage 10 summaries and Stage 07 metadata plus an injected renderer stub. Assert creation of `primary.haplotype_network_summary.tsv`, `primary.haplotype_network_commands.tsv`, and `primary.haplotype_network_report.md`. The report must contain both organelles, `pegas::haploNet`, complete-case filtering, and “not ancestry proportions.”

- [ ] **Step 2: Verify RED**

Expected: `run_haplotype_network_analysis` is absent.

- [ ] **Step 3: Implement the stage API**

Add frozen `HaplotypeNetworkInput` and `HaplotypeNetworkResult` dataclasses and:

```python
def run_haplotype_network_analysis(
    snp_alignment_dir=DEFAULT_SNP_ALIGNMENT_DIR,
    metadata_path=DEFAULT_METADATA_PATH,
    output_dir=DEFAULT_OUTPUT_DIR,
    run_label="primary",
    rscript=Path("Rscript"),
    renderer_path=DEFAULT_RENDERER_PATH,
    runner=subprocess.run,
):
    inputs = read_haplotype_network_inputs(snp_alignment_dir, run_label)
    metadata = read_sample_metadata(metadata_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    commands = []
    for network_input in inputs:
        records = read_fasta(network_input.alignment_fasta_path)
        validate_sample_metadata([sample_id for sample_id, _ in records], metadata)
        filtered, kept, dropped = filter_complete_case_sites(records)
        paths = network_paths(output_dir, network_input.organelle, run_label)
        write_network_input_fasta(paths.input_fasta, filtered)
        write_network_site_table(paths.site_table, read_tsv(network_input.site_table_path), kept)
        write_network_metadata(paths.metadata_table, filtered, metadata)
        paths.popart_nexus.write_text(build_popart_nexus(filtered, metadata))
        command = build_renderer_command(rscript, renderer_path, paths.input_fasta, paths.metadata_table, paths.prefix, network_input.organelle)
        completed = runner(command, text=True, capture_output=True, check=False)
        if completed.returncode:
            raise HaplotypeNetworkError(completed.stderr or completed.stdout)
        validate_renderer_outputs(paths)
        result = read_renderer_result(network_input, paths, len(kept), len(dropped))
        results.append(result)
        commands.append(command)
    write_haplotype_network_outputs(output_dir, results, commands, run_label)
    return results
```

The implementation reads the Stage 10 summary, validates sample identity/order, prepares each organelle, writes PopART NEXUS, runs R, validates products, reads renderer counts, writes combined summary/command/report files, and returns two results. Store exact commands with `shlex.join`.

- [ ] **Step 4: Add the runner wrapper**

```python
#!/usr/bin/env python3
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from dudleya_organelle_alignment_pipeline.haplotype_network import main
if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 5: Run targeted and full tests, then commit**

Run the Stage 21 tests and full discovery suite with the pipeline environment. Commit with message `feat: add Stage 21 haplotype network analysis`.

## Task 6: Generate and inspect primary results

**Files:**
- Create: `dudleya_organelle_alignment_pipeline/results/21_haplotype_network/*`

- [ ] **Step 1: Generate results**

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" R_DEFAULT_DEVICE=png python3 dudleya_organelle_alignment_pipeline/scripts/run_haplotype_network.py --run-label primary
```

Expected counts: cpDNA 275 samples, 1,977 retained sites, 151 haplotypes; mtDNA 275 samples, 116 retained sites, 58 haplotypes.

- [ ] **Step 2: Validate all products**

Assert 275 unique assignments per organelle; assignment frequencies sum to 275; retained+dropped sites equal 2,015 and 146; all edge endpoints exist; every required table and figure is nonempty.

- [ ] **Step 3: Inspect original-resolution PNGs**

Check titles, legends, node sizes, consistent species colors, mutation labels, clipping, and readability. If cpDNA is crowded, adjust canvas size and label policy only; do not change the network method.

- [ ] **Step 4: Commit results**

Commit `results/21_haplotype_network` with message `results: add primary organelle haplotype networks`.

## Task 7: Document Stage 21 and final deliverables

**Files:**
- Modify: `dudleya_organelle_alignment_pipeline/PROCESS.md`
- Modify: `dudleya_organelle_alignment_pipeline/README.md`
- Modify: `dudleya_organelle_alignment_pipeline/results/organelle_population_report.md`
- Modify: `dudleya_organelle_alignment_pipeline/results/final_deliverables_manifest.tsv`
- Regenerate: `dudleya_organelle_alignment_pipeline/results/13_tool_audit/*`

- [ ] **Step 1: Update PROCESS and README**

Change the canonical range to 00-21. Document inputs, complete-case filtering, `pegas::haploNet`, the runner command, outputs, species-sector meaning, PopART NEXUS export, and limitations. Retain Stage 16/18 documentation unchanged.

- [ ] **Step 2: Update the integrated report**

Add “Haploid Haplotype Networks” after PCA. Populate actual counts and links. State that nodes/edges are not populations, ancestry proportions, or known ancestral transitions. Retain the complete ADMIXTURE section.

- [ ] **Step 3: Extend the deliverables manifest**

Append rows for both network PNGs, both assignment TSVs, both PopART NEXUS files, and the combined Stage 21 report. Do not remove or edit existing rows.

- [ ] **Step 4: Regenerate and check the tool audit**

Run `run_tool_audit.py --audit-label primary`; confirm `r_pegas` is `FOUND` and total tool count rises by one.

- [ ] **Step 5: Commit documentation**

Commit documentation, manifest, and regenerated audit with message `docs: document organelle haplotype network deliverables`.

## Task 8: Final verification

- [ ] **Step 1: Run all tests fresh**

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" python3 -m unittest discover -s dudleya_organelle_alignment_pipeline/tests -v
```

Expected: all existing 70 tests plus all new Stage 21 tests pass.

- [ ] **Step 2: Re-run Stage 21**

Run the primary Stage 21 command again and confirm deterministic counts and nonempty outputs.

- [ ] **Step 3: Verify scope**

Use `git status --short` and the implementation commit file lists. Expected: clean status and no changed path under `dudleya_organelle_alignment_pipeline/study_guide/`.

- [ ] **Step 4: Verify deliverables manually**

Confirm every new manifest path exists, all PNG/PDF/SVG files open, the Stage 21 report names `pegas` and complete-case filtering, and the integrated report contains both ADMIXTURE and haplotype-network sections.
