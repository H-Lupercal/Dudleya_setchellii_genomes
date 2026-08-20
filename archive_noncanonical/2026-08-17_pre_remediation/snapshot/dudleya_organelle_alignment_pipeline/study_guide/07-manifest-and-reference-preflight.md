# Chapter 7 — Manifest and Reference Preflight (Stages 00–01)

> Part 2 of 4 · Pipeline Walkthrough · Prev:
> [Organelle Biology](./06-organelle-biology.md) · Next: [Pilot Mapping and
> Investigations](./08-pilot-mapping-and-investigations.md)

Every pipeline chapter follows the same eight steps: the question, the files, the
runner-to-module map, the code, the Python concepts, the biological result, the
failure modes, and exercises. This first one goes slowly because it also
establishes the reading rhythm.

## 7.1 The question

Before a single read is aligned, you must answer boring questions correctly, or
everything downstream inherits the mistake: *Which FASTQ files exist? Which R1
and R2 files belong to the same biological sample? Which sample uses which naming
convention? Which samples have population metadata? Which are safe to analyze,
and which must be excluded?* Then: *Is the reference structurally what we expect,
are the tools installed, and which small, diverse pilot set should we test on
first?*

Stages 00 and 01 do only this. They **do not** align reads, trim, call variants,
or build consensus. Getting identity and pairing right up front is what makes the
expensive stages trustworthy. `[CODE]`

## 7.2 The files

**Stage 00** ([`manifest.py`](../manifest.py), run by
[`../scripts/build_sample_manifest.py`](../scripts/build_sample_manifest.py)):

- Input: FASTQ files under `genomicsDrive_data_dump/`, and the population-code
  CSV `genomicsDrive_data_dump/QB3.Berkeley.251217/Dudleya DNAx - Population
  Codes.csv`.
- Output in `results/00_manifest/`: `samples.tsv` (one row per sample),
  `analysis_samples.tsv` (the primary paired-end set), `excluded_samples.tsv`,
  `pairing_report.tsv`, and
  [`preflight_summary.md`](../results/00_manifest/preflight_summary.md).

**Stage 01** ([`prepare_reference_and_pilot.py`](../prepare_reference_and_pilot.py),
run by
[`../scripts/prepare_reference_and_pilot.py`](../scripts/prepare_reference_and_pilot.py)):

- Input: the combined reference
  `../../dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa`
  and Stage 00's `analysis_samples.tsv`.
- Output in `results/01_reference_pilot/`: `reference_checks.tsv`,
  `tool_checks.tsv`, `index_checks.tsv`, `pilot_samples.tsv`, and
  `reference_pilot_summary.md`. It also creates the `.fai` and `bwa index` files
  next to the reference, but only if the tools are installed.

## 7.3 Runner → module → external commands → results

Stage 00 shells out to nothing — it is pure Python filesystem work. Stage 01
optionally runs `samtools faidx` and `bwa index`
([Chapter 4](./04-shell-and-external-tools.md)), but only after checking the
tools are present. Both stages' runner scripts are the standard thin wrappers
([Chapter 3, §3.1](./03-reusable-code-patterns.md)).

## 7.4 The code: parsing a FASTQ filename into identity

The heart of Stage 00 is turning a filename into structured metadata. The parser
recognizes three naming profiles that share the Illumina suffix
`_S<n>_L<lane>_R<1|2>_001.fastq.gz`. The regular expression that peels off that
suffix is:

```python
FASTQ_NAME_RE = re.compile(
    r"^(?P<prefix>.+)_"
    r"(?P<sequencing_sample>S\d+)_"
    r"(?P<lane>L\d+)_"
    r"(?P<read>R[12])_"
    r"(?P<chunk>\d+)"
    r"(?P<extension>\.(?:fastq|fq)(?:-\d+)?(?:\.gz)?)$",
    re.IGNORECASE,
)
```

A **regular expression** ("regex") is a pattern-matching mini-language.
`(?P<name>...)` captures a named group you can pull out later; `\d+` means "one
or more digits"; `R[12]` means "R followed by 1 or 2". Everything before the
`_S...` suffix is captured as `prefix`, which is the biological sample ID. The
`(?:-\d+)?` allows chunked names like `..._001.fastq-011.gz`, and a test
confirms discovery still finds those. `[TEST]`

The prefix is then classified into one of three profiles by `classify_prefix`:

```python
def classify_prefix(prefix: str) -> tuple[str, str, str, str]:
    main_match = MAIN_STANDARD_RE.match(prefix)
    if main_match:
        return ("main_standard", main_match.group("popcode"),
                main_match.group("du_id"), main_match.group("lp_id"))
    du_lp_match = INITIAL_DU_LP_RE.match(prefix)
    if du_lp_match:
        return ("initial_du_lp", "", du_lp_match.group("du_id"),
                du_lp_match.group("lp_id"))
    du_dash_match = INITIAL_DU_DASH_RE.match(prefix)
    if du_dash_match:
        return ("initial_du_dash", "", du_dash_match.group("du_id"), "")
    return ("unrecognized", "", "", "")
```

- **`main_standard`** — e.g. `CY_RED_LP_202_Du-561`: has a population code
  (`CY_RED`), a plant ID (`LP_202`), and a DNA ID (`Du-561`). Only these carry
  population metadata.
- **`initial_du_dash`** — e.g. `DU-4A`: an early batch with no encoded population.
- **`initial_du_lp`** — e.g. `DU014LP012`: another early batch, DU and LP but no
  population code.

Tests parse one example of each profile and assert every field, so the regexes
are pinned to real filenames — including a tricky variant where LP and Du are
separated by a hyphen instead of an underscore (`CY_ALA_LP_298-Du-767`). `[TEST]`

### From records to pairs

`build_manifest` groups records by `(batch, sample_id)`, splits them into R1 and
R2, and decides a `pair_status`:

```python
def determine_pair_status(r1_records, r2_records) -> str:
    if not r1_records and not r2_records:
        return "missing_R1_and_R2"
    if not r1_records:
        return "missing_R1"
    if not r2_records:
        return "missing_R2"
    if len(r1_records) != len(r2_records):
        return "uneven_read_counts"
    if len(r1_records) > 1:
        return "complete_multi_file"
    return "complete"
```

Only `complete` (exactly one R1 and one R2) becomes
`include_primary_paired_end`; missing mates become `exclude_missing_mate`;
anything else is `review_before_primary_analysis`. That translation is
`determine_analysis_status`, and it also writes the human explanation string that
travels with the sample all the way to the exclusion audit. A test builds a
three-file fixture (one complete pair, one R1-only sample) and asserts the
complete sample is included while the lone R1 is flagged `missing_R2` and
excluded. `[TEST]`

### Population metadata

`load_population_codes` reads the population CSV into a code-keyed dict, tolerating
the real-world messiness that the code column is titled
`Code (if it doesn't start with a TWO letter code = DUSE)`:

```python
code_field = next(
    (field for field in reader.fieldnames if field.lower().startswith("code")),
    None,
)
```

`next((... for ...), None)` returns the first column whose name starts with
"code", or `None`. When a `main_standard` sample's popcode is in the CSV, its
species and population name are filled and `metadata_status` becomes `resolved`;
otherwise the species is *inferred* from the popcode prefix
(`infer_species_from_popcode`: `CY_` → *D. cymosa*, `AB` → *D. abramsii*,
otherwise *D. setchellii*). The initial batches stay
`unresolved_initial_sample`. `[CODE]`

## 7.5 The Python concepts here

- **Regex with named groups** for filename parsing (new in this chapter).
- **`defaultdict(list)`** to bucket FASTQ records by sample without pre-checking
  keys.
- **`try/except ValueError`** around `parse_fastq_path`: an unparseable name is
  caught and recorded as a `ManifestIssue`, not allowed to crash the run — the
  pipeline records the problem and continues.
- **`write_dataclass_tsv`** deriving columns from the `ManifestRow` dataclass
  ([Chapter 3, §3.2](./03-reusable-code-patterns.md)).
- **`next(iterable, default)`** to find the first match safely.

## 7.6 Stage 01: validating the reference and choosing a pilot

`validate_reference_records` confirms the reference contains exactly the expected
records with the expected lengths:

```python
EXPECTED_REFERENCE_LENGTHS = {"chloroplast": 150274, "mitochondria": 243359}
```

A missing or extra record raises `ReferenceValidationError`; a length mismatch is
a non-fatal `WARN`. A test proves a missing `mitochondria` record raises. `[TEST]`
`check_tools` records whether `bwa`, `samtools` (required) and `fastp`, `fastqc`,
`multiqc`, `bcftools` (recommended) are on `PATH`, and `prepare_indexes` runs
`samtools faidx` and `bwa index` **only** for the tools that are present — so the
stage degrades gracefully on a machine without the aligner installed rather than
crashing. `[CODE]`

`select_pilot_samples` then chooses up to 15 diverse, complete samples
deterministically. It first guarantees one representative of each early batch and
each of the three species, then round-robins across species to add more
populations, always skipping excluded samples and de-duplicating population
codes. The observed pilot is 15 samples: 2 initial-batch representatives, 5
*D. cymosa*, 3 *D. abramsii*, 5 *D. setchellii*. `[RESULT]` Three tests check the
selection logic: it skips excluded rows, prefers main-dataset populations after
the initial representatives, and round-robins across species groups. `[TEST]`

Why deterministic? Because a pilot you cannot reproduce is not a pilot. The
selection uses sorted keys throughout (`sample_sort_key`, `population_sort_key`),
so the same inputs always yield the same 15 samples.

## 7.7 The biological result, and what it does *not* claim

`preflight_summary.md` reports counts by batch, naming profile, pair status,
metadata status, and analysis status. The honest reading:

- **What it establishes:** which files exist, which pair, which have metadata,
  and which 278 complete pairs are eligible for the primary analysis. `[RESULT]`
- **What it does not:** it says nothing about sequencing *quality* or organelle
  *coverage* — a sample can be "complete" here and still be dropped later for low
  coverage. Pairing is necessary, not sufficient. `[CODE]`

The two manually verified missing-mate samples
(`ABAB_MAD_LP_225_Du-592` missing R1, `QUI1_LP_256_Du-655` missing R2) are kept
in `samples.tsv` for the audit trail but excluded from `analysis_samples.tsv`.
The policy is explicit: if such a sample is ever aligned single-end, that must be
a separate sensitivity analysis, never mixed into the primary paired-end
dataset. `[CODE]`

## 7.8 Failure modes

- **Unparseable FASTQ name** → recorded as `unparsed_fastq_name` in
  `pairing_report.tsv`; the run continues. `[CODE]`
- **Missing mate** → `missing_R1`/`missing_R2`, excluded from the analysis set
  with a note. `[CODE]`
- **Uneven or multi-file reads** → `review_before_primary_analysis`; the pilot
  and alignment stages will *raise* if they later meet a multi-file row they were
  not told to expect ([Chapter 8](./08-pilot-mapping-and-investigations.md)).
- **Missing population CSV** → `load_population_codes` returns an empty dict;
  metadata is simply unresolved rather than crashing. `[CODE]`
- **Wrong reference records** → `ReferenceValidationError`, halting Stage 01
  before any indexing. `[CODE]`
- **Missing aligner** → indexing is *skipped* with a recorded reason, not a
  crash; the missing tool is listed in the summary's "Before Pilot Alignment"
  section. `[CODE]`

## 7.9 Exercises

1. **Trace.** Given `ABMU_HOR_LP_140_Du-410_S77_L006_R2_001.fastq.gz`, what are
   `sample_id`, `read`, `naming_profile`, and `popcode` after
   `parse_fastq_path`? Which species would `infer_species_from_popcode` assign if
   the popcode is absent from the CSV?
2. **Predict.** A sample has two R1 files and two R2 files (a re-sequenced
   library). What `pair_status` does `determine_pair_status` return, and what
   `analysis_status`? Will the pilot stage accept it?
3. **Predict.** You delete the population CSV and rerun Stage 00. What happens to
   the `metadata_status` of the `CY_RED_LP_202_Du-561` sample, and does the run
   crash?
4. **Modify.** Suppose a new batch uses names like `DUSE_412_S9_L001_R1_001...`.
   Which regex would you add, and where in `classify_prefix` would you slot it so
   the existing profiles still match first?
5. **Debug.** A colleague reports that `analysis_samples.tsv` is missing a sample
   they know is complete. List three checks — filename, batch grouping, pair
   status — you would run, and which output file confirms each.
6. **Interpret.** The preflight says 278 complete pairs, but the final analysis
   uses 275 samples. Is that a bug in Stage 00? Explain what removed the other
   three and why Stage 00 is not the place that would catch them.

Solutions are in [Chapter 19](./19-solutions.md).

> Next: [Chapter 8 — Pilot Mapping and Reference Investigations (Stages 02–04)](./08-pilot-mapping-and-investigations.md)
