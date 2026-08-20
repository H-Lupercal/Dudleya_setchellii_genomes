# Chapter 14 — The Bioinformatics Tool Audit (Stage 13)

> Part 2 of 4 · Pipeline Walkthrough · Prev:
> [Population Fst](./13-population-fst.md) · Next: [Reading the
> Trees](./15-phylogenetics-interpretation.md)

Stage 13 answers a plain operational question — *is the software this pipeline
needs actually installed?* — and it is the cleanest example in the repo of a
software-engineering pattern worth learning: **dependency injection for
testability.**

## 14.1 The question

*Which required and recommended tools are present, at what versions, and is the
machine ready to reproduce the completed pipeline and to run the remaining
analyses?*

It depends on nothing upstream — you can run it any time — which is why it sits
off to the side of the data-flow map ([Chapter 1](./01-data-flow-map.md)).

## 14.2 The files

[`tool_audit.py`](../tool_audit.py) (runner:
[`../scripts/run_tool_audit.py`](../scripts/run_tool_audit.py)) writes
`results/13_tool_audit/primary.tool_audit.tsv` and
`primary.tool_audit_report.md`. The TSV is the machine-readable record; the
Markdown report is its human-readable summary and is what the top-level README
points to.

## 14.3 The code: specs, necessity, and how a "tool" is defined

The list of tools is data, not logic — a tuple of `ToolSpec` records:

```python
@dataclass(frozen=True)
class ToolSpec:
    tool_id: str
    executables: tuple[str, ...]        # accept any of these names
    necessity: str
    required_for: str
    version_args: tuple[str, ...]        # how to ask for a version
```

Necessity has three levels that map to the project's phases:

- **`required_current`** — needed to reproduce the completed mapping/QC/variant/
  tree pipeline (`python3`, `bwa`, `samtools`, `bcftools`, `fastp`, `fastqc`,
  `multiqc`, `iqtree`).
- **`required_remaining`** — needed for the remaining analyses (`plink`,
  `admixture`, `Rscript`, and the Python libraries `matplotlib`, `pandas`,
  `sklearn`, `biopython`, plus R `ggplot2`, `ape`).
- **`recommended_remaining`** — nice to have (`FastTree`, `vcftools`, `bedtools`,
  `seaborn`, `ete3`, `patchwork`, `snakemake`).

The clever part is that a "tool" can be a **Python library checked by importing
it**. The spec for matplotlib runs `python3 -c "import matplotlib;
print(matplotlib.__version__)"`:

```python
ToolSpec("python_matplotlib", ("python3",), "required_remaining",
         "PCA scatterplots, tree rendering, and static figures",
         ("-c", "import matplotlib; print(matplotlib.__version__)")),
```

So "is matplotlib installed?" is answered by actually importing it, not by
finding a binary — exactly right, because a bare `python3` on `PATH` says nothing
about whether the plotting stages will work. A test confirms all the
visualization dependencies are in the default specs. `[TEST]`

## 14.4 The pattern to learn: dependency injection

`check_tool` does the work, and its signature is the lesson:

```python
def check_tool(spec, resolver=shutil.which, runner=run_version_command) -> ToolResult:
    for executable in spec.executables:
        path = resolver(executable)
        if path:
            try:
                version = first_version_line(runner([path, *spec.version_args])) ...
            except Exception as exc:
                return ToolResult(..., status="MISSING", ...,
                                  note=f"Executable found, but version/import check failed: {exc}")
            return ToolResult(..., status="FOUND", ...)
    return ToolResult(..., status="MISSING", note=missing_note(spec.necessity))
```

`resolver` (how to find a tool) and `runner` (how to run it) are **parameters
with real defaults**: in production they are `shutil.which` and the real
subprocess runner. But a test can pass *fake* ones and check the logic without
any tools installed:

```python
result = check_tool(spec,
    resolver=lambda executable: "/env/bin/samtools",
    runner=lambda command: "samtools 1.23.1\nUsing htslib")
self.assertEqual(result.status, "FOUND")
self.assertEqual(result.version, "samtools 1.23.1")
```

That is **dependency injection**: pass the risky, environment-dependent behavior
in as an argument so it can be swapped for a stub in tests. It is why the audit's
tests run on any machine and still cover the found path, the missing path, and
the subtle "executable exists but the import fails" path (a `python3` that lacks
matplotlib is reported MISSING with a clear note). `[TEST]` If you learn one
reusable engineering idea from this pipeline, make it this one.

## 14.5 Summarizing readiness

`summarize_audit` counts what is missing by necessity, and two `@property`
methods turn that into go/no-go answers:

```python
@property
def ready_for_current_pipeline(self) -> bool:
    return not self.missing_required_current

@property
def ready_for_remaining_goal(self) -> bool:
    return not self.missing_required_current and not self.missing_required_remaining
```

The report prints an explicit interpretation line for each — "the completed
pipeline can be reproduced" or "do not continue until missing tools are
installed." A test builds a two-tool audit (iqtree present, admixture missing)
and asserts `ready_for_remaining_goal` is `False` with `admixture` in the missing
list. `[TEST]`

## 14.6 The result

The project's audit records the pinned versions from [Chapter 4,
§4.1](./04-shell-and-external-tools.md) as FOUND — `bwa` 0.7.19, `samtools`/
`bcftools` 1.23, `IQ-TREE` 3.1.2, `PLINK` 1.9, `ADMIXTURE` 1.3.0, and the Python
and R plotting stacks — so the environment is ready for both the completed and
remaining analyses. `[RESULT]` This is a reproducibility artifact: it lets anyone
check *before* rerunning whether their machine can produce the same outputs.

## 14.7 The Python concepts here

- **Data-driven design**: the tool list is a tuple of dataclasses, so adding a
  tool is adding a row, not writing code.
- **Dependency injection** via default-valued function parameters.
- **`Callable` type hints** (`Resolver = Callable[[str], str | None]`) documenting
  the injected function shapes.
- **`@property` booleans** turning counts into decisions.
- **Broad `except Exception`** *on purpose* here, to convert any import/version
  failure into a recorded MISSING rather than a crash — the rare justified use of
  a wide catch.

## 14.8 Failure modes

- **A required tool missing** → recorded MISSING with an "Install before..." note;
  the report's readiness line flips to not-ready. No crash — reporting the gap is
  the job. `[CODE]`
- **Executable present but broken** (e.g. `python3` without `sklearn`) → MISSING
  with the failing message captured. `[CODE]`
- **Version flag differs across tool versions** → `first_version_line` keeps just
  the first non-empty line, tolerating chatty output. `[CODE]`

## 14.9 Exercises

1. **Trace.** `check_tool` is given a spec for `plink` with `executables =
   ("plink", "plink2")`. If `resolver` returns `None` for `plink` but a path for
   `plink2`, what happens?
2. **Predict.** A `python_sklearn` spec runs on a machine with `python3` but no
   scikit-learn. What `status` and `note` does `check_tool` return?
3. **Predict.** All `required_current` tools are found but `admixture` is missing.
   What do `ready_for_current_pipeline` and `ready_for_remaining_goal` return?
4. **Modify.** You want to audit `mafft`. Write the `ToolSpec` line, choosing a
   necessity level and version args, and say where it goes in `TOOL_SPECS`.
5. **Debug.** The audit reports `iqtree` MISSING on a machine where `iqtree2`
   runs fine. What is the fix — the spec, the resolver, or the `PATH`?
6. **Interpret.** Why does the pipeline check libraries by importing them instead
   of trusting that `python3` on `PATH` is enough? Give a concrete stage that
   would break if it trusted `python3` alone.

Solutions in [Chapter 19](./19-solutions.md).

> Part 2 complete. Next: [Chapter 15 — Reading the Phylogenetic Trees](./15-phylogenetics-interpretation.md)
