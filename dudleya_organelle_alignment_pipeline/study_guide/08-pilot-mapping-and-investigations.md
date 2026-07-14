# Chapter 8 — Pilot Mapping and Reference Investigations (Stages 02–04)

> Part 2 of 4 · Pipeline Walkthrough · Prev:
> [Manifest and Reference Preflight](./07-manifest-and-reference-preflight.md) ·
> Next: [Masks, Alignment, and Sample QC](./09-masks-alignment-and-sample-qc.md)

## 8.1 The question

Two questions, really. First: *does organelle mapping even work on these data,
and how strong is the cpDNA and mtDNA signal?* You answer that on a small,
diverse pilot before committing to all 275 samples. Second, and more subtle:
*which parts of each organelle can we actually trust?* The pilot mapping surfaces
two problems — mitochondrial repeats and the chloroplast inverted repeat — that
Stages 03 and 04 investigate and turn into evidence for the masks in Stage 05.

Stage 02 aligns the 15 pilot samples and summarizes mapping. It does **not** call
variants, build consensus, make final alignments, or run any population analysis.
`[CODE]`

## 8.2 The files

**Stage 02** ([`pilot_alignment.py`](../pilot_alignment.py), run by
[`../scripts/run_pilot_alignment.py`](../scripts/run_pilot_alignment.py)):

- Input: `results/01_reference_pilot/pilot_samples.tsv` and the combined
  reference (with its indexes from Stage 01).
- Output in `results/02_pilot_alignment/`:
  `pilot_alignment_sample_summary.tsv` (one row per sample),
  `pilot_alignment_by_organelle.tsv` (one row per sample × organelle),
  `pilot_alignment_report.md`, `commands.tsv`, and the generated
  `bam/`, `qc/`, `logs/` (git-ignored analysis artifacts).

**Stages 03–04** are *investigations*, not scripted modules. Their evidence lives
in `results/03_mtdna_investigation/` and `results/04_cpdna_investigation/` and is
consumed directly by Stage 05. The key files are
`mtdna_high_mapq_consensus_intervals.tsv` and `cpdna_self_repeat_intervals.tsv`.

## 8.3 Runner → module → external commands → results

Stage 02 runs the `bwa mem | samtools view | samtools sort` pipe plus the four
QC commands per sample, exactly as taught in [Chapter 4,
§4.2](./04-shell-and-external-tools.md). It records each command in
`commands.tsv` and each tool's stderr in a per-sample log. This is where the
alignment functions the whole pipeline reuses are *defined* — Stage 06 imports
them wholesale ([Chapter 3, §3.7](./03-reusable-code-patterns.md)).

## 8.4 The code: from BAM to organelle metrics

After a sample's BAM exists, Stage 02 parses two `samtools` outputs.
`parse_idxstats_file` reads mapped-read counts per reference record, skipping the
`*` unmapped row:

```python
def parse_idxstats_file(path: Path) -> dict[str, int]:
    ...
    if len(fields) < 4 or fields[0] == "*":
        continue
    mapped_counts[fields[0]] = int(fields[2])   # column 3 = mapped reads
```

So a BAM becomes `{"chloroplast": 42, "mitochondria": 7}` — the per-organelle
read split, from one shared BAM. A test asserts exactly this on a hand-written
idxstats fixture. `[TEST]`

`parse_depth_file` scans the per-base depth file and accumulates, per organelle,
the total depth and the number of bases at ≥1×, ≥5×, ≥10×:

```python
depth = int(fields[2])
counters[organelle]["total_depth"] += depth
if depth >= 1:  counters[organelle]["bases_ge_1x"]  += 1
if depth >= 5:  counters[organelle]["bases_ge_5x"]  += 1
if depth >= 10: counters[organelle]["bases_ge_10x"] += 1
```

Those raw counts become an `OrganelleMetrics` dataclass whose `@property`
methods derive `mean_depth` and `breadth_ge_Nx` on demand ([Chapter 2,
§2.3](./02-python-essentials.md)). Because `samtools depth -aa` emits every
position including zeros, the breadth denominator is the full reference length —
"breadth at 1×" is literally the fraction of the reference covered by at least
one read. A test checks a tiny depth file yields mean 2.5 and breadth 0.5 for a
4-base record. `[TEST]`

### The QC flags

`build_sample_summary` turns the metrics into a human-readable row and, crucially,
a list of QC notes triggered by two thresholds:

```python
LOW_MAPPED_READS_THRESHOLD = 100
LOW_BREADTH_THRESHOLD = 0.50
...
if cp_mapped < LOW_MAPPED_READS_THRESHOLD:    notes.append("low_chloroplast_mapped_reads")
if mt_mapped < LOW_MAPPED_READS_THRESHOLD:    notes.append("low_mitochondria_mapped_reads")
if cp_metrics.breadth_ge_1x < LOW_BREADTH_THRESHOLD: notes.append("low_chloroplast_breadth_ge_1x")
if mt_metrics.breadth_ge_1x < LOW_BREADTH_THRESHOLD: notes.append("low_mitochondria_breadth_ge_1x")
```

A sample with no flags gets the note `pass_initial_mapping_screen`. A test feeds
a sample with 1000 cp reads and 2 mt reads and asserts
`low_mitochondria_mapped_reads` fires while the cp fraction is 0.998004. `[TEST]`
These are *screening* thresholds, meant to catch obviously broken samples early;
the stricter downstream QC is Stage 06.

## 8.5 The Python concepts here

- **Streaming subprocess pipes** with `Popen` and `.wait()`, checking every exit
  code ([Chapter 3, §3.5](./03-reusable-code-patterns.md)).
- **Atomic output**: write to `sample.tmp.bam`, then `tmp_bam.replace(bam_path)`
  so a partial BAM never looks finished.
- **`@property`-derived metrics** keeping raw counts and computed rates separate.
- **Resumability** via `outputs_are_ready` + `--force`/`--refresh-qc`
  ([Chapter 3, §3.6](./03-reusable-code-patterns.md)).
- **A hand-rolled `median_breadth`** that sorts values and averages the middle
  two for even counts — a reminder that the pipeline avoids heavy dependencies
  for simple stats.

## 8.6 The result, and the two problems it exposed

The pilot summarized 15 samples, 30 sample-by-organelle rows, tens of millions of
mapped reads, a median input organelle mapping fraction around 0.144, median
chloroplast breadth ≥1× of ~0.9999, and median mitochondrial breadth ≥1× of
~0.96. `[RESULT]` Read that carefully: chloroplast is essentially fully covered,
mitochondria nearly so *at permissive mapping quality*. But that mtDNA number is
a trap, and finding the trap is the point of Stage 03.

### Stage 03 — mtDNA repeats

The report notes that "after correcting the `samtools depth` quality flags"
([Chapter 4, §4.3](./04-shell-and-external-tools.md)), the remaining mtDNA issue
is repeat/ambiguity handling: a **high-MAPQ-only** check shows much lower
*unique-placement* breadth than the permissive-MAPQ depth. In plain terms, many
mtDNA positions look covered only because reads that could belong to several
repeat copies were counted; restrict to confidently, uniquely placed reads and
much of the mitochondrion falls away. The investigation writes
`mtdna_high_mapq_consensus_intervals.tsv` — the intervals that survive at high
MAPQ across enough pilot samples — which becomes the mtDNA population track in
Stage 05. `[RESULT]` `[BIO]`

### Stage 04 — the chloroplast inverted repeat

The cpDNA verification confirms the expected large inverted repeat and writes
`cpdna_self_repeat_intervals.tsv` recording the two IR copies (normalized
coordinates ~`82091–107826` and ~`124539–150274`). The conclusion is that
all-sample chloroplast processing is safe *if* you keep only one IR copy for
population genetics. `[RESULT]` This evidence becomes
`cpdna_population_sites.bed` in Stage 05.

This is the pipeline being honest with itself: the pilot did not just say "yes,
mapping works." It said "mapping works, *and here are the two regions you must
not trust naively*," and it turned that judgment into machine-readable interval
tables rather than prose.

## 8.7 Failure modes

- **Missing tools or reference indexes** → `require_tools` /
  `require_reference_indexes` raise before any alignment, naming what to install
  or re-run. `[CODE]`
- **A row marked `complete` contains multiple R1 or R2 paths** →
  `read_alignment_samples` raises `AlignmentError`, because this stage expects
  exactly one R1 and one R2. Rows whose `pair_status` is not `complete` are
  skipped before this check. `[TEST]`
- **A truncated FASTQ** → `count_fastq_records` raises on a line count not
  divisible by four. `[CODE]`
- **A subprocess in the pipe fails** → any non-zero return among bwa/view/sort
  raises `AlignmentError` naming the log to inspect. `[CODE]`
- **Biological failure: permissive-MAPQ over-optimism.** The subtle one — mtDNA
  looks well-covered until you demand unique placement. This is not a crash; it
  is a scientific trap the investigation exists to catch, and ignoring it would
  produce false mtDNA variants downstream. `[BIO]`

## 8.8 Exercises

1. **Trace.** A pilot sample has idxstats rows `chloroplast 150274 5000 0`,
   `mitochondria 243359 40 0`, `* 0 0 900`. What does `parse_idxstats_file`
   return, and which QC note(s) will `build_sample_summary` add?
2. **Predict.** You rerun Stage 02 with `--refresh-qc`. Which files are
   regenerated and which are reused? What `commands.tsv` step names appear?
3. **Predict.** Two positions have depth 4 and 6 in a 2-base reference. What are
   `mean_depth`, `breadth_ge_1x`, and `breadth_ge_5x`?
4. **Modify.** You want to add a ≥20× breadth metric. Which fields of
   `OrganelleMetrics` and which lines of `parse_depth_file` and
   `build_organelle_summary_rows` change?
5. **Debug.** A sample shows `chloroplast_breadth_ge_1x = 0.999` but
   `mitochondria_breadth_ge_1x = 0.02`. Is this a mapping failure, a biological
   fact, or expected? What would you check in the mtDNA investigation before
   deciding?
6. **Interpret.** The pilot reports median mtDNA breadth ≥1× of 0.96, yet Stage
   05 keeps only ~44,930 bp of the 243,359 bp mitochondrion for variants.
   Reconcile these two facts using the permissive-vs-high-MAPQ distinction.

Solutions in [Chapter 19](./19-solutions.md).

> Next: [Chapter 9 — Masks, All-Sample Alignment, and Sample QC (Stages 05–07)](./09-masks-alignment-and-sample-qc.md)
