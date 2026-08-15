# Concatenated Organelle Bootstrap Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create and verify a 275-sample cpDNA-then-mtDNA callable-consensus alignment, while preparing but not executing a 10,000-replicate combined IQ-TREE runner.

**Architecture:** Add an isolated concatenation module and CLI that validate paired FASTA records, concatenate by sample identifier in cpDNA order, and emit both provenance metadata and the minimal summary schema consumed by the unchanged phylogenetic-tree engine. Copy the existing fixed 10,000-replicate wrapper to a combined-analysis wrapper with isolated Stage 22 input and Stage 23 output paths. Existing scripts and results remain untouched.

**Tech Stack:** Python 3 standard library, `unittest`, existing `dudleya_organelle_alignment_pipeline.phylogenetic_tree` engine, IQ-TREE 3.1.2 (prepared only; not run in this plan).

---

## File Structure

- Create `dudleya_organelle_alignment_pipeline/concatenated_consensus.py`: FASTA validation, per-sample concatenation, output metadata, CLI argument handling.
- Create `dudleya_organelle_alignment_pipeline/scripts/run_concatenated_consensus.py`: thin executable wrapper around the new module.
- Create `dudleya_organelle_alignment_pipeline/scripts/run_concatenated_phylogenetic_tree_10000.py`: copied fixed runner pointing at combined input/output directories.
- Create `dudleya_organelle_alignment_pipeline/tests/test_concatenated_consensus.py`: exact seam, ordering, validation, and output tests.
- Create `dudleya_organelle_alignment_pipeline/tests/test_run_concatenated_phylogenetic_tree_10000.py`: fixed bootstrap-runner configuration test.
- Create at runtime `dudleya_organelle_alignment_pipeline/results/22_concatenated_consensus/`: verified concatenated FASTA and metadata only.

### Task 1: Concatenation validation and data model

**Files:**
- Create: `dudleya_organelle_alignment_pipeline/tests/test_concatenated_consensus.py`
- Create: `dudleya_organelle_alignment_pipeline/concatenated_consensus.py`

- [ ] **Step 1: Write failing tests for exact concatenation and cpDNA ordering**

Add tests that write cpDNA records `S1=ACGN`, `S2=TTAA` and reverse-ordered mtDNA records `S2=GG`, `S1=NC`, then call:

```python
alignment = concatenate_consensus_alignments(cp_path, mt_path)
assert alignment.sample_names == ("S1", "S2")
assert alignment.sequences == {"S1": "ACGNNC", "S2": "TTAAGG"}
assert alignment.cpdna_length == 4
assert alignment.mtdna_length == 2
assert alignment.combined_length == 6
assert alignment.mtdna_start == 5
```

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```bash
python3 -m unittest dudleya_organelle_alignment_pipeline.tests.test_concatenated_consensus -v
```

Expected: import failure because `concatenated_consensus.py` does not exist.

- [ ] **Step 3: Implement the minimal validated data model**

Implement:

```python
class ConcatenatedConsensusError(RuntimeError): ...

@dataclass(frozen=True)
class FastaAlignment:
    sample_names: tuple[str, ...]
    sequences: dict[str, str]
    sequence_length: int

@dataclass(frozen=True)
class ConcatenatedAlignment:
    sample_names: tuple[str, ...]
    sequences: dict[str, str]
    cpdna_length: int
    mtdna_length: int

    @property
    def combined_length(self) -> int:
        return self.cpdna_length + self.mtdna_length

    @property
    def mtdna_start(self) -> int:
        return self.cpdna_length + 1

    @property
    def missing_bases(self) -> int:
        return sum(sequence.count("N") for sequence in self.sequences.values())
```

Implement `read_fasta_alignment(path)` to preserve header order, uppercase sequence text, reject sequence text before a header, empty alignments, duplicate identifiers, empty sequences, and inconsistent record lengths. Implement `concatenate_consensus_alignments(cpdna_path, mtdna_path)` to reject unequal identifier sets and construct each record as the exact cpDNA sequence plus the mtDNA sequence with the same identifier.

- [ ] **Step 4: Run the focused test and verify GREEN**

Run the Task 1 unittest command. Expected: exact-concatenation and ordering tests pass.

- [ ] **Step 5: Add validation tests one behavior at a time**

Add individual tests asserting `ConcatenatedConsensusError` for:

```python
duplicate FASTA identifiers
mismatched cpDNA/mtDNA identifier sets
inconsistent cpDNA record lengths
inconsistent mtDNA record lengths
empty FASTA input
```

Run each newly added test first to verify it fails for the intended reason, then add the smallest validation needed and rerun until it passes.

- [ ] **Step 6: Run all concatenation model tests**

Run the Task 1 unittest command. Expected: all Task 1 tests pass with no warnings.

- [ ] **Step 7: Commit Task 1**

```bash
git add dudleya_organelle_alignment_pipeline/concatenated_consensus.py dudleya_organelle_alignment_pipeline/tests/test_concatenated_consensus.py
git commit -m "Add validated organelle consensus concatenation"
```

### Task 2: Concatenation outputs and CLI

**Files:**
- Modify: `dudleya_organelle_alignment_pipeline/concatenated_consensus.py`
- Modify: `dudleya_organelle_alignment_pipeline/tests/test_concatenated_consensus.py`
- Create: `dudleya_organelle_alignment_pipeline/scripts/run_concatenated_consensus.py`

- [ ] **Step 1: Write failing output tests**

Add a temporary-directory test calling:

```python
result = run_concatenation(
    cpdna_path=cp_path,
    mtdna_path=mt_path,
    output_dir=output_dir,
    run_label="primary",
)
```

Assert that:

```python
result.sample_count == 2
result.combined_length == 6
read_fasta_alignment(result.fasta_path).sequences["S1"] == "ACGNNC"
```

Assert `primary.callable_consensus_summary.tsv` contains exactly one combined row with fields `organelle=cpDNA_mtDNA`, `track_id=cpdna_then_mtdna`, `sample_count=2`, `consensus_length=6`, `missing_bases`, and `alignment_fasta_path`. Assert the detailed summary and report contain `cpDNA_end=4` and `mtDNA_start=5`.

- [ ] **Step 2: Run the output test and verify RED**

Run the Task 1 unittest command. Expected: failure because `run_concatenation` and output writers are absent.

- [ ] **Step 3: Implement deterministic output writing**

Add `ConcatenationResult`, `write_fasta`, `write_tsv`, `write_outputs`, and `run_concatenation`. Use defaults:

```python
DEFAULT_CPDNA_PATH = Path("dudleya_organelle_alignment_pipeline/results/11_callable_consensus/cpDNA.primary.callable_consensus.fa")
DEFAULT_MTDNA_PATH = Path("dudleya_organelle_alignment_pipeline/results/11_callable_consensus/mtDNA.primary.callable_consensus.fa")
DEFAULT_OUTPUT_DIR = Path("dudleya_organelle_alignment_pipeline/results/22_concatenated_consensus")
DEFAULT_RUN_LABEL = "primary"
```

Write FASTA records in cpDNA order with 80-character sequence lines. Write the compatibility summary with only the six fields required by `phylogenetic_tree.read_tree_inputs`. Write a detailed summary containing input paths, sample count, component lengths, combined length, 1-based boundaries, missing-base counts, and output path. Write a Markdown report stating that mtDNA was appended unchanged to cpDNA.

- [ ] **Step 4: Add the thin CLI wrapper**

Copy the import/bootstrap structure of `scripts/run_phylogenetic_tree.py`, importing `main` from `concatenated_consensus`. Make it executable and accept optional input/output arguments through the module parser.

- [ ] **Step 5: Run tests and CLI help**

```bash
python3 -m unittest dudleya_organelle_alignment_pipeline.tests.test_concatenated_consensus -v
python3 dudleya_organelle_alignment_pipeline/scripts/run_concatenated_consensus.py --help
```

Expected: all tests pass; help lists `--cpdna-path`, `--mtdna-path`, `--output-dir`, and `--run-label`.

- [ ] **Step 6: Commit Task 2**

```bash
git add dudleya_organelle_alignment_pipeline/concatenated_consensus.py dudleya_organelle_alignment_pipeline/scripts/run_concatenated_consensus.py dudleya_organelle_alignment_pipeline/tests/test_concatenated_consensus.py
git commit -m "Add concatenated consensus outputs and runner"
```

### Task 3: Prepared combined 10,000-bootstrap runner

**Files:**
- Create: `dudleya_organelle_alignment_pipeline/scripts/run_concatenated_phylogenetic_tree_10000.py`
- Create: `dudleya_organelle_alignment_pipeline/tests/test_run_concatenated_phylogenetic_tree_10000.py`

- [ ] **Step 1: Write the failing fixed-arguments test**

Load the runner with `importlib.util.spec_from_file_location`, call `build_run_arguments()`, and assert:

```python
value("--consensus-dir") == "dudleya_organelle_alignment_pipeline/results/22_concatenated_consensus"
value("--output-dir") == "dudleya_organelle_alignment_pipeline/results/23_concatenated_bootstrap_phylogenetic_tree_10000"
value("--bootstrap-replicates") == "10000"
value("--threads") == "14"
value("--run-label") == "primary"
```

- [ ] **Step 2: Run the runner test and verify RED**

```bash
python3 -m unittest dudleya_organelle_alignment_pipeline.tests.test_run_concatenated_phylogenetic_tree_10000 -v
```

Expected: file-not-found failure because the copied runner does not exist.

- [ ] **Step 3: Copy and specialize the existing runner**

Copy `scripts/run_phylogenetic_tree_10000.py` to `scripts/run_concatenated_phylogenetic_tree_10000.py`. Preserve the original. Change only the docstring, fixed input/output constants, and returned arguments. Continue calling the unchanged `phylogenetic_tree.main`.

- [ ] **Step 4: Run the runner and phylogenetic command tests**

```bash
python3 -m unittest \
  dudleya_organelle_alignment_pipeline.tests.test_run_concatenated_phylogenetic_tree_10000 \
  dudleya_organelle_alignment_pipeline.tests.test_phylogenetic_tree -v
```

Expected: all tests pass. Do not execute the runner itself.

- [ ] **Step 5: Commit Task 3**

```bash
git add dudleya_organelle_alignment_pipeline/scripts/run_concatenated_phylogenetic_tree_10000.py dudleya_organelle_alignment_pipeline/tests/test_run_concatenated_phylogenetic_tree_10000.py
git commit -m "Prepare concatenated 10000-bootstrap runner"
```

### Task 4: Full verification and controlled concatenation

**Files:**
- Runtime outputs: `dudleya_organelle_alignment_pipeline/results/22_concatenated_consensus/`

- [ ] **Step 1: Run all focused tests fresh**

```bash
python3 -m unittest \
  dudleya_organelle_alignment_pipeline.tests.test_concatenated_consensus \
  dudleya_organelle_alignment_pipeline.tests.test_run_concatenated_phylogenetic_tree_10000 \
  dudleya_organelle_alignment_pipeline.tests.test_phylogenetic_tree -v
```

Expected: zero failures and zero errors.

- [ ] **Step 2: Run only the concatenation CLI**

```bash
python3 dudleya_organelle_alignment_pipeline/scripts/run_concatenated_consensus.py
```

Expected: reports 275 samples, cpDNA length 124,538, mtDNA length 44,930, and combined length 169,468. Do not invoke IQ-TREE.

- [ ] **Step 3: Verify all 275 seams against source data**

Use the project Python environment and Biopython to load the three FASTAs. Assert:

```python
len(combined) == len(cpDNA) == len(mtDNA) == 275
set(combined) == set(cpDNA) == set(mtDNA)
combined[sample] == cpDNA[sample] + mtDNA[sample]  # every sample
len(combined[sample]) == 169468                    # every sample
```

Expected: print `verified_exact_seams=275`, `combined_length=169468`, and no assertion failures.

- [ ] **Step 4: Verify metadata and absence of a bootstrap process**

Confirm detailed/compatibility summaries contain one combined row and the exact boundary. Check process state for `run_concatenated_phylogenetic_tree_10000` and IQ-TREE; expected: no matching bootstrap process started by this task.

- [ ] **Step 5: Check repository state and hand off**

Run `git diff --check` and `git status --short`. Report the concatenation output paths, exact 275/275 seam verification, combined length, test count, and prepared runner path. Explicitly state that the 10,000-replicate bootstrap has not started and awaits user approval.
