# Chapter 3 — Reusable Code Patterns

> Part 1 of 4 · Foundations · Prev: [Python
> Essentials](./02-python-essentials.md) · Next: [Shell and External
> Tools](./04-shell-and-external-tools.md)

The pipeline has sixteen implementation modules but only about seven patterns. Learn them once
here and every stage chapter gets shorter, because it can say "this is the
standard runner wrapper" or "this is the standard resumable subprocess call"
instead of re-explaining. When you open a new module, spend your attention on
what is *different*; the scaffolding is almost always one of these.

## 3.1 The runner-script wrapper (`scripts/run_*.py`)

Every stage has a tiny script in [`../scripts/`](../scripts/). They are all the
same shape. Here is the entire `scripts/build_sample_manifest.py`:

```python
#!/usr/bin/env python3
"""Run: build the Dudleya organelle FASTQ manifest."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from dudleya_organelle_alignment_pipeline.manifest import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
```

That is the *whole* file, and every runner is identical except for which
module's `main` it imports. What it does:

- `Path(__file__).resolve().parents[2]` walks up two directories from the script
  (`scripts/` → pipeline dir → repository root) and calls that `ROOT`.
- `sys.path.insert(0, str(ROOT))` puts the repository root on Python's import
  path, so `import dudleya_organelle_alignment_pipeline...` works no matter what
  directory you run from.
- It imports the module's `main` and runs it.

**Why separate the runner from the module at all?** So the logic lives in an
importable module the tests can call directly, while the script is a thin,
path-fixing entry point for humans. The `# noqa: E402` comment tells the linter
"yes, this import is deliberately below the `sys.path` line." Once you have seen
one runner, you have seen all sixteen. `[CODE]`

## 3.2 `read_tsv` and `write_tsv`

You already met `read_tsv` in [Chapter 2, §2.7](./02-python-essentials.md). Its
partner writes a list of dicts back out with a fixed column order:

```python
def write_tsv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
```

You pass an explicit `fieldnames` list so the columns come out in a stable,
documented order rather than whatever order the dict happened to have. Almost
every module carries its own copy of this pair — deliberately, so each stage is
a self-contained unit you can read without chasing imports.

There is a dataclass-aware variant in the earliest modules,
`write_dataclass_tsv`, which derives the column names from the dataclass fields
themselves:

```python
def write_dataclass_tsv(path: Path, rows: list[object], row_type: type[object]) -> None:
    fieldnames = list(row_type.__dataclass_fields__.keys())
    ...
```

`row_type.__dataclass_fields__` is the automatically generated map of a
dataclass's fields; taking its keys gives the column order for free. This is how
[`manifest.py`](../manifest.py) writes `samples.tsv` straight from
`ManifestRow` objects. `[CODE]`

## 3.3 The "read summary → build inputs → run → write outputs" shape

From Stage 08 onward, every stage has the same four-beat rhythm, because each
stage's input is the *previous stage's summary table*:

1. **Read the upstream summary TSV** to learn where the real data files are.
2. **Validate** that those files exist (and raise if not).
3. **Do the work** (call a tool, or compute in pure Python).
4. **Write** this stage's own summary TSV, a `report.md`, and the data product.

For example, [`variant_filtering.py`](../variant_filtering.py) reads
`primary.variant_calling_summary.tsv` to find the raw VCFs, filters them, and
writes `primary.variant_filtering_summary.tsv`, which
[`snp_alignment.py`](../snp_alignment.py) then reads to find the filtered VCFs.
The summary TSVs are the *contracts between stages*. When you want to know what
a stage consumes, open the `read_*_inputs` function at the top of its module; it
always names the summary file and the columns it trusts. `[CODE]`

This is why you rarely pass data in memory between stages: each stage is a
separate process, run separately, that communicates only through files on disk.
That makes every stage independently rerunnable and independently auditable.

## 3.4 Run labels: `labeled_output_name`

Filenames like `primary.variant_calling_summary.tsv` get their prefix from one
small function defined in [`variant_calling.py`](../variant_calling.py) and
imported by everything downstream:

```python
def labeled_output_name(name: str, run_label: str) -> str:
    if not run_label:
        return name
    return f"{run_label}.{name}"
```

So with `run_label="primary"`, `variant_calling_summary.tsv` becomes
`primary.variant_calling_summary.tsv`; with an empty label it stays unprefixed.
There is a sibling, `label_output_prefix`, that inserts the label *in the middle*
of a name so `cpDNA.raw` becomes `cpDNA.primary.raw` (used for the per-organelle
VCF prefixes). Between them, every downstream file carries the run label, which
is exactly what keeps a `smoke` test run from ever colliding with the `primary`
deliverables. `[CODE]` This is the mechanism behind the run-label discussion in
[Chapter 1, §1.4](./01-data-flow-map.md).

## 3.5 Driving external tools with `subprocess`

The pipeline runs command-line tools through Python's `subprocess` module in two
styles. The simple style, for a single command whose output goes to a file:

```python
def run_command(args: list[str]) -> None:
    subprocess.run(args, check=True)
```

`subprocess.run(args, check=True)` runs a command given as a **list of strings**
(never one big string — that avoids shell-quoting bugs) and, with `check=True`,
raises if the command exits non-zero.

The advanced style builds a **pipeline of processes**, mirroring a shell's
`bwa mem ... | samtools view ... | samtools sort ...`. From
[`pilot_alignment.py`](../pilot_alignment.py):

```python
bwa_proc = subprocess.Popen(align_command, stdout=subprocess.PIPE, stderr=log_handle)
view_proc = subprocess.Popen(view_command, stdin=bwa_proc.stdout,
                             stdout=subprocess.PIPE, stderr=log_handle)
bwa_proc.stdout.close()
sort_proc = subprocess.Popen(sort_command, stdin=view_proc.stdout,
                             stdout=subprocess.DEVNULL, stderr=log_handle)
view_proc.stdout.close()
sort_return = sort_proc.wait()
view_return = view_proc.wait()
bwa_return = bwa_proc.wait()
```

`Popen` starts a process without waiting. Wiring one process's `stdout` into the
next's `stdin` builds the pipe entirely in memory, so the aligned reads are
never written to a temporary SAM file — they stream straight through the filter
and the sort. Each `.wait()` returns the exit code; the code checks all three so
a failure anywhere in the pipe is caught. `stderr=log_handle` captures every
tool's diagnostics into a per-sample log. The details of what these particular
commands *do* are in [Chapter 4](./04-shell-and-external-tools.md). `[CODE]`

The mapping, variant-calling, filtering, and tree stages record their external
commands verbatim in a `commands.tsv` (via the `shlex_join` helper, which safely
quotes arguments). ADMIXTURE records its invocations in the ADMIXTURE summary
table instead. These provenance records are first-class deliverables, not debug
artifacts.

## 3.6 Resumability: `outputs_are_ready`, `--force`, `--refresh-qc`

The heavy stages (mapping, variant calling) can take hours across 275 samples,
so they are **resumable**: if a sample's outputs already exist, the stage reuses
them instead of recomputing. The pattern is a small `outputs_are_ready`-style
check plus `--force`:

```python
if force or not outputs.bam_path.exists():
    # (re)run the alignment
elif refresh_qc or not outputs_are_ready(outputs):
    # keep the BAM, but regenerate index/flagstat/idxstats/depth
else:
    # record "reuse_existing_outputs" and move on
```

So a crashed run can be restarted and it will skip everything already finished.
`--force` recomputes from scratch; `--refresh-qc` keeps the expensive alignment
but rebuilds the cheap QC files (useful when only the depth-quality flags
changed). Variant calling and filtering use the same idea keyed on the output
VCF's existence and non-zero size. When you see a `"reuse_existing_outputs"` row
in a `commands.tsv`, that is this pattern reporting a skip. `[CODE]`

## 3.7 Modules that build on other modules

The pipeline does not copy the alignment machinery twice. Stage 06
([`all_sample_alignment.py`](../all_sample_alignment.py)) *imports* the core
functions from Stage 02 ([`pilot_alignment.py`](../pilot_alignment.py)):

```python
from dudleya_organelle_alignment_pipeline.pilot_alignment import (
    AlignmentError, AlignmentSample, build_sample_summary, count_fastq_records,
    outputs_for_sample, parse_depth_file, parse_idxstats_file,
    read_alignment_samples, run_alignment_commands, run_qc_commands, ...
)
```

The pilot stage is where the alignment functions are *defined and first tested*;
the all-sample stage *reuses* them and adds only the track-aware coverage
summary that is new at that scale. Likewise, several late stages import
`read_fasta` and `read_sample_metadata` from
[`pca_analysis.py`](../pca_analysis.py), and everything downstream of variant
calling imports `labeled_output_name` from
[`variant_calling.py`](../variant_calling.py). When a chapter says "the same
function you met in Stage 02," this is why — it is literally the same function.
`[CODE]`

A quick map of the important cross-module reuse:

| Function | Defined in | Reused by |
|---|---|---|
| `safe_sample_name`, `shlex_join`, alignment/QC helpers | `pilot_alignment.py` | `all_sample_alignment.py`, `variant_calling.py`, `callable_consensus.py` |
| `labeled_output_name` | `variant_calling.py` | filtering, SNP, consensus, trees, PCA, admixture, popgen, tree viz |
| `count_vcf_records` | `variant_calling.py` | `variant_filtering.py` |
| `read_fasta`, `read_sample_metadata`, `choose_plot_group` | `pca_analysis.py` | `admixture_analysis.py`, `population_genetics.py` |
| `EXPECTED_REFERENCE_LENGTHS` | `prepare_reference_and_pilot.py` | `analysis_masks.py` |

## 3.8 Tests as executable specifications

Each module has a matching test file in [`../tests/`](../tests/), and the tests
follow their own repeated shapes: build a tiny fixture in a
`tempfile.TemporaryDirectory()`, call one function, and assert on the exact
result. Two kinds dominate:

- **Command-builder tests** assert the *exact argument list* a function produces
  without running the tool — e.g. `test_build_iqtree_command_...` checks that
  `--fast` is present for the initial run and `-B 1000 --bnni` for the final
  run. These pin the tool invocations down to the flag. `[TEST]`
- **Pure-logic tests** feed in a hand-written table, VCF, or FASTA and assert the
  computed output byte for byte — e.g. the callable-consensus test asserts the
  literal strings `"CTTNNG"` and `"CGTANG"` come out of a six-base region
  ([Chapter 23](./23-capstone-sample-trace.md) traces that one by hand).

Because the tests never invoke the heavy bioinformatics binaries (`bwa`,
`iqtree`, `admixture`, `plink`), the suite runs in about a second. There are
70 tests across the 16 files. Most need only the standard library, but the PCA
rendering test imports NumPy, scikit-learn, and matplotlib, while two tree-rendering
tests import Biopython and matplotlib. Run the complete suite in the pipeline
environment so those Python dependencies are available. The tests are the
fastest, truest description of what each function is supposed to do. `[TEST]`

```bash
# Run all 70 with the pipeline environment:
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 -m unittest discover -s dudleya_organelle_alignment_pipeline/tests -v
```

## 3.9 The pattern checklist

When you open any module, you can now expect:

- [ ] A module docstring stating what the stage does and, often, what it
  deliberately does *not* do.
- [ ] `DEFAULT_*` constants naming inputs, outputs, and thresholds.
- [ ] Frozen dataclasses for the records and results.
- [ ] A custom `*Error` exception class.
- [ ] `read_*_inputs` reading an upstream summary TSV.
- [ ] The core work in one or two functions.
- [ ] `write_*_outputs` producing a summary TSV plus a `report.md`.
- [ ] `build_arg_parser` + `main` + the `if __name__ == "__main__"` line.

Everything else is the biology and the specific computation, which is what the
stage chapters are for.

> Next: [Chapter 4 — The Shell and External Tools](./04-shell-and-external-tools.md)
