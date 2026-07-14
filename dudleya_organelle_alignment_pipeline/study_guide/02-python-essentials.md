# Chapter 2 — Python Essentials for This Codebase

> Part 1 of 4 · Foundations · Prev: [Data-Flow
> Map](./01-data-flow-map.md) · Next: [Reusable Code
> Patterns](./03-reusable-code-patterns.md)

This chapter teaches only the Python you need to *read* this pipeline. It is not
a general Python course. Every example is a real construct pulled from the
repository, so when you meet it again in a stage chapter, it will already be
familiar. If you know these ten things, no line in the pipeline will stop you.

The declared environment targets modern Python (3.11+) and leans on a small, consistent set
of features. Let us take them in the order you meet them when you open a module.

## 2.1 `from __future__ import annotations`

Every implementation module starts with this line:

```python
from __future__ import annotations
```

It changes one thing: type hints are treated as text, not evaluated at import
time. The practical effect for a reader is that you can write a type like
`int | str` or `list[Path]` even on older interpreters, and you can reference a
class in its own method signatures before it is fully defined. You do not need
to do anything with this line; just know it makes the type hints below legal and
free. `[CODE]`

## 2.2 Type hints

Type hints annotate what a function expects and returns. They do not enforce
anything at runtime — Python will not stop you passing the wrong type — but they
document intent precisely, which is exactly what you want when reading unfamiliar
code. From [`manifest.py`](../manifest.py):

```python
def determine_pair_status(
    r1_records: list[FastqRecord],
    r2_records: list[FastqRecord],
) -> str:
```

Read this as: "takes two lists of `FastqRecord` objects, returns a string." A
few notations you will see constantly:

- `list[Path]` — a list of `Path` objects.
- `dict[str, int]` — a dict mapping strings to ints.
- `tuple[int, int]` — a 2-element tuple of ints (a coordinate pair).
- `int | None` — "an int, or `None`." The `| None` marks an optional value.
- `set[str] | None = None` — an optional set that defaults to `None`.

When a signature says `-> None`, the function is called for its side effects
(writing a file, running a command), not for a return value.

## 2.3 Dataclasses: the pipeline's data records

A **dataclass** is a class whose main job is to hold named fields. Instead of
writing an `__init__` by hand, you list the fields with their types and Python
generates the constructor for you. The pipeline uses dataclasses everywhere to
model rows and results. From [`manifest.py`](../manifest.py):

```python
@dataclass(frozen=True)
class FastqRecord:
    path: Path
    filename: str
    batch: str
    sample_id: str
    naming_profile: str
    sequencing_sample: str
    lane: str
    read: str
    chunk: str
    popcode: str = ""
    du_id: str = ""
    lp_id: str = ""
```

Three things to notice:

1. **`@dataclass`** is a *decorator* — a line starting with `@` that modifies the
   class below it. Here it generates `__init__`, so you can write
   `FastqRecord(path=..., filename=..., ...)`.
2. **`frozen=True`** makes instances immutable: once created, you cannot reassign
   a field. This is a deliberate safety choice — a parsed FASTQ record should
   never be quietly mutated later. Attempting `record.read = "R2"` raises an
   error. `[CODE]`
3. **Defaults** like `popcode: str = ""` mean those fields are optional at
   construction. Fields with defaults must come after fields without them.

You access fields with a dot: `record.sample_id`, `record.read`. That is all a
dataclass is — a labeled container with a free constructor.

### Computed fields with `@property`

Some dataclasses expose values that are *computed* from their fields rather than
stored. From [`pilot_alignment.py`](../pilot_alignment.py):

```python
@dataclass(frozen=True)
class OrganelleMetrics:
    organelle: str
    reference_length: int
    total_depth: int
    bases_ge_1x: int
    bases_ge_5x: int
    bases_ge_10x: int

    @property
    def mean_depth(self) -> float:
        if self.reference_length == 0:
            return 0.0
        return self.total_depth / self.reference_length

    @property
    def breadth_ge_1x(self) -> float:
        return self._breadth(self.bases_ge_1x)
```

A `@property` is a method you call *without parentheses*, as if it were a field:
you write `metrics.mean_depth`, not `metrics.mean_depth()`. It recomputes on
every access. The pipeline uses this so a metrics object stores only the raw
counts (`total_depth`, `bases_ge_1x`) and derives means and breadths on demand,
guarding against divide-by-zero. `[CODE]` The `_breadth` method has a leading
underscore, a convention meaning "internal helper, not part of the public
surface."

## 2.4 `pathlib.Path`: filenames as objects

The pipeline never manipulates paths as raw strings. It uses `Path` objects from
the `pathlib` module. A `Path` knows how to join, inspect, read, and write
itself:

```python
from pathlib import Path

bam_path = output_dir / "bam" / f"{stem}.organelle.sorted.bam"   # join with /
bam_path.exists()               # bool: does the file exist?
bam_path.parent                 # the directory containing it
bam_path.name                   # just the filename
bam_path.stem                   # filename without the last suffix
bam_path.suffix                 # ".bam"
bam_path.with_suffix(".tmp.bam")# swap the suffix
bam_path.as_posix()             # forward-slash string, for command lines
path.mkdir(parents=True, exist_ok=True)   # create dir (and parents), no error if present
path.write_text("...")          # write a whole string to the file
path.read_text()                # read the whole file as a string
```

The `output_dir / "bam" / filename` syntax reads oddly at first — the `/`
operator joins path segments. It is cleaner and more portable than string
concatenation, and it is used on nearly every line that touches the filesystem.
`.as_posix()` shows up whenever a path is about to become an argument to an
external command, because command-line tools want a plain string. `[CODE]`

## 2.5 f-strings

An **f-string** is a string prefixed with `f` where `{...}` is replaced by the
value inside. The pipeline builds nearly every filename, message, and report
line this way:

```python
marker_id = f"{admixture_input.organelle}_snp_{site_index + 1}"
# -> "cpDNA_snp_1"

f"[{index}/{len(samples)}] {sample.sample_id}"
# -> "[3/15] CY_RED_LP_202_Du-561"
```

You can format numbers inside the braces. `f"{value:.6f}"` prints six decimal
places; `f"{variance[0] * 100:.2f}%"` prints a percentage with two decimals.
This is how the reports get their tidy `0.999993` and `36.62%` figures. `[CODE]`

## 2.6 Comprehensions and generators

A **list comprehension** builds a list in one expression. A **generator
expression** does the same lazily, one item at a time, without building the
whole list in memory. You will read dozens of both. From
[`variant_calling.py`](../variant_calling.py):

```python
# list comprehension: keep only rows whose track is wanted
tracks = {row["track_id"]: row for row in read_tsv(track_table)
          if row["track_id"] in wanted}          # this is a *dict* comprehension

# generator inside a function call: sum without an intermediate list
total_mapped = sum(int(row["total_organelle_mapped_reads"])
                   for row in sample_summaries)
```

The shapes to recognize:

- `[expr for item in iterable]` — a list.
- `{key: value for item in iterable}` — a dict.
- `{expr for item in iterable}` — a set (deduplicates).
- `(expr for item in iterable)` — a generator; common as the sole argument to
  `sum(...)`, `sorted(...)`, `any(...)`, `all(...)`, or `"".join(...)`.

Add `if condition` at the end to filter. This one idiom — build a filtered
collection in a single readable line — is the backbone of every "read a table,
keep the rows I want" step in the pipeline.

## 2.7 Reading and writing TSV with `csv`

Every table in this pipeline is a tab-separated file (TSV). The `csv` module
reads and writes them. The pipeline uses `DictReader` (each row becomes a dict
keyed by column name) and `DictWriter` (each dict becomes a row). This helper
appears, nearly identically, in almost every module:

```python
import csv

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
```

`with path.open(...) as handle:` opens the file and guarantees it is closed
afterward, even if an error is raised — that is what `with` does. `DictReader`
turns the header row into dict keys, so downstream code reads `row["sample_id"]`
instead of counting columns. `delimiter="\t"` says tab, not comma. Every value
comes back as a **string**, which is why you see `int(row["sample_count"])` and
`float(row["fst"])` conversions all over the pipeline. `[CODE]` The mirror-image
writer is covered in [Chapter 3, §3.2](./03-reusable-code-patterns.md).

## 2.8 Dicts, sets, `defaultdict`, and `Counter`

Beyond plain dicts and sets, the pipeline uses two specialized containers from
`collections`:

- **`defaultdict(list)`** — a dict that auto-creates an empty list the first time
  you touch a missing key, so you can write `groups[popcode].append(...)`
  without checking whether `popcode` is present yet. Used to bucket samples by
  population in [`population_genetics.py`](../population_genetics.py).
- **`Counter`** — a dict that counts things. `Counter("ACGT")` counts letters;
  `counts1 + counts2` merges two counters by summing. It is the heart of the
  allele-frequency math in Fst ([Chapter 13](./13-population-fst.md)). `[CODE]`

Sets matter for a subtle reason you will see in the mask code: `{a, b} - {b}`
computes a set difference. The failed-variant mask in
[`callable_consensus.py`](../callable_consensus.py) is literally
`raw_site_keys - filtered_site_keys` — "sites present in the raw calls but not in
the filtered calls." Set algebra makes that a one-liner.

## 2.9 Exceptions and custom error classes

When a stage cannot safely continue, it *raises* an exception rather than
returning a bad result. Each module defines its own exception class so the
failure is self-describing:

```python
class VariantCallingError(RuntimeError):
    """Raised when this stage cannot safely call variants."""

# ...later...
if not bam_path.exists():
    raise VariantCallingError(f"Missing BAM for {sample_id}: {bam_path}")
```

`class VariantCallingError(RuntimeError)` means "a new kind of error that *is a*
`RuntimeError`." Raising it stops the stage immediately with a message naming
the exact file. You will see `AlignmentError`, `MaskDefinitionError`,
`CallableConsensusError`, `PcaAnalysisError`, and siblings — one per module. The
tests deliberately trigger these (`with self.assertRaises(AlignmentError):`) to
prove the guardrails fire ([Chapter 3, §3.6](./03-reusable-code-patterns.md)).
This "fail loudly and early" style is central to a pipeline that must be *right*,
not just finish. `[CODE]`

`try/except` appears where the pipeline can *recover* or *translate* an error —
for example, catching an unparseable FASTQ name and recording it as an issue
instead of crashing ([`manifest.py`](../manifest.py) `build_manifest`), or
catching a missing optional library and re-raising it as a clear
"install this tool" message ([`pca_analysis.py`](../pca_analysis.py) `run_pca`).

## 2.10 `argparse` and the `main(argv)` convention

Every runnable module ends with an `argparse`-based command-line interface and a
`main` function:

```python
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="...")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--sample-id", action="append", dest="sample_ids")
    return parser

def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    ...
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
```

Reading argparse definitions tells you a stage's knobs at a glance:

- `type=int` converts the text argument to an int; `default=4` is used if the
  flag is absent.
- `action="store_true"` makes a boolean flag: `--force` present means `True`.
- `action="append"` lets a flag repeat: `--sample-id A --sample-id B` collects
  `["A", "B"]`. `dest="sample_ids"` renames where it lands on `args`.
- A flag written `--min-mapq` becomes `args.min_mapq` (dash to underscore).

`main(argv=None)` takes its arguments as a parameter so tests can call it
directly, and `argparse` falls back to the real command line when `argv` is
`None`. `if __name__ == "__main__":` means "only run this when the file is
executed directly, not when it is imported." `raise SystemExit(main())` turns
`main`'s integer return into the process exit code (0 = success). `[CODE]`

## 2.11 One more: `zip(..., strict=True)`

`zip` pairs up two sequences element by element. The pipeline pairs VCF sample
columns with sample names this way, and passes `strict=True` so that a
length mismatch *raises* instead of silently truncating:

```python
for sample, genotype_field in zip(sample_names, columns[9:], strict=True):
    sequence_parts[sample].append(genotype_to_base(genotype_field, ref, alt))
```

If the VCF ever had a different number of sample columns than sample names, that
would be a serious corruption — `strict=True` guarantees you hear about it rather
than getting a quietly wrong alignment. `[CODE]` This is the same defensive
instinct as the custom exceptions: prefer a loud failure to a plausible-looking
wrong answer.

## 2.12 What you can now read

With dataclasses, `@property`, `Path`, f-strings, comprehensions, `csv`
dict-based I/O, the `collections` helpers, custom exceptions, argparse, and
`zip(strict=True)`, you have the entire Python vocabulary of this repository. The
only remaining building blocks are *patterns* built from these primitives —
runner wrappers, the TSV writer, run-label naming, subprocess piping, and
resumability — which the next chapter covers once so the stage chapters do not
have to.

> Next: [Chapter 3 — Reusable Code Patterns](./03-reusable-code-patterns.md)
