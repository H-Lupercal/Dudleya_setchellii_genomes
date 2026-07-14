# Chapter 9 — Masks, All-Sample Alignment, and Sample QC (Stages 05–07)

> Part 2 of 4 · Pipeline Walkthrough · Prev: [Pilot
> Mapping and Investigations](./08-pilot-mapping-and-investigations.md) · Next:
> [Variants to Alignments](./10-variants-to-alignments.md)

These three stages turn pilot judgment into documented analysis tracks (05), scale
mapping to every sample with track-aware QC (06), and freeze the exact 275-sample
set the rest of the pipeline uses (07). This is where "QC tracks are not
population tracks" stops being advice and becomes BED files and code.

## 9.1 The question

*Which regions of each organelle do we use, for what purpose? Do all 275+ samples
sequence well enough to keep? What is the final, frozen sample set?*

## 9.2 Stage 05 — analysis masks

[`analysis_masks.py`](../analysis_masks.py) (runner:
[`../scripts/build_analysis_masks.py`](../scripts/build_analysis_masks.py))
converts the Stage 03/04 evidence into six BED tracks plus two manifest tables.
It does **not** align or call variants — it only defines regions. `[CODE]`

Inputs: `cpdna_self_repeat_intervals.tsv` (Stage 04) and
`mtdna_high_mapq_consensus_intervals.tsv` (Stage 03). Outputs in
`results/05_analysis_masks/`: six `.bed` files, `analysis_regions.tsv` (the
interval audit in both coordinate systems), `analysis_tracks.tsv` (the
machine-readable contract), and `mask_summary.md`.

### The six tracks

| Track | Organelle | Purpose | Notes |
|---|---|---|---|
| `cpdna_full_coverage` | cpDNA | sample QC | whole 150,274 bp reference |
| `cpdna_ir_regions` | cpDNA | annotation | both IR copies, documented |
| `cpdna_duplicate_ir_mask` | cpDNA | mask | the later IR copy, to be excluded |
| `cpdna_population_sites` | cpDNA | variants + popgen | full reference minus the duplicate IR copy → **124,538 bp** |
| `mtdna_permissive_coverage` | mtDNA | sample QC | whole 243,359 bp reference |
| `mtdna_high_confidence_unique` | mtDNA | variants + popgen | high-MAPQ intervals → **44,930 bp** |

`[RESULT]` The design in one sentence: **coverage-QC tracks span the whole
organelle; population-genetic tracks span only the regions you can trust.**

### The code: keeping one IR copy

`build_cpdna_tracks` reads the two largest self-repeat intervals (the IR copies),
marks the later copy as the duplicate mask, and builds the population track as
the *complement* of that mask:

```python
population_regions=complement_regions(
    track_id="cpdna_population_sites", ...,
    reference_length=reference_length,
    masked_regions=duplicate_mask,   # the later IR copy
    ...)
```

`complement_regions` walks the reference from position 1, emitting every stretch
*not* covered by a mask. `read_major_cpdna_repeat_pair` requires exactly two
major IR intervals and raises `MaskDefinitionError` otherwise — a guard against
malformed repeat evidence. A test builds a tiny 20-bp reference with IR copies at
`5–8` and `15–18`, then asserts the population sites come out as `1–14` and
`19–20` (everything except the masked second copy). `[TEST]` That is the
inverted-repeat biology of [Chapter 6, §6.3](./06-organelle-biology.md) made
concrete.

### The code: the mtDNA threshold

`build_mtdna_tracks` keeps only the high-MAPQ intervals supported by at least
`DEFAULT_MTDNA_HIGH_CONFIDENCE_THRESHOLD = 12` usable pilot samples:

```python
rows = [row for row in read_tsv(high_mapq_intervals_path)
        if int(row["threshold_usable_samples"]) == threshold]
```

Requiring an interval to be uniquely mappable in a dozen pilot samples is what
shrinks the trusted mtDNA to ~44,930 bp. A test with a `12`-sample interval and a
`10`-sample interval confirms only the `12`-sample intervals are kept. `[TEST]`

### The contract: `analysis_tracks.tsv`

This TSV is the authoritative machine-readable statement of each track's intended
purpose and downstream use. Stage 06 carries `purpose` into its QC summaries but
does not reject tracks by purpose. Stage 08 selects the two expected population
tracks by their hard-coded `track_id` values and does not inspect `purpose`.
Consequently, the purpose column is documentation rather than a fully enforced
type system; editing track identities or meanings can still create a biologically
wrong run. A test asserts the six tracks appear in the expected order. `[TEST]`

## 9.3 Stage 06 — all-sample alignment with track-aware QC

[`all_sample_alignment.py`](../all_sample_alignment.py) (runner:
[`../scripts/run_all_sample_alignment.py`](../scripts/run_all_sample_alignment.py))
maps every primary paired-end sample from `analysis_samples.tsv` and summarizes
coverage three ways: per organelle, and per analysis track.

It imports the alignment and QC machinery from Stage 02 unchanged ([Chapter 3,
§3.7](./03-reusable-code-patterns.md)) and adds only the *track* summary. The new
piece is `parse_track_depth_file`, which counts depth **only inside** the BED
intervals of each track:

```python
for region in regions_by_record.get(record, []):
    if region.start_1based <= position <= region.end_1based:
        counter = counters[region.track_id]
        counter["total_depth"] += depth
        if depth >= 1: counter["bases_ge_1x"] += 1
        ...
```

The denominator for a track's breadth is the total BED-defined length, so a
position never observed in the depth file counts as zero coverage — exactly right
for asking "what fraction of this trusted region did the sample cover?" Tests
feed a depth file spanning inside and outside the track and assert only in-track
positions are counted, and that breadth comes out as `0.666667` for a 3-bp region
with 2 covered bases. `[TEST]`

### The result and the QC decision

Stage 06 writes per-sample, per-organelle, and per-track summaries, plus
`downstream_sample_qc_decisions.tsv`. Reviewing that QC, the pipeline excludes
three samples — `ABAB_MAD_LP_222_Du-589`, `CY_HUN_LP_265_Du-684`,
`CY_RED_LP_202_Du-561` — the three lowest-input samples, which failed one or both
organelle coverage screens. `[RESULT]` The command is resumable: an existing BAM
and its QC are reused unless `--force` or `--refresh-qc` is passed, and the report
is rewritten after every sample so a long run always has an up-to-date summary.

## 9.4 Stage 07 — freezing the downstream set

[`downstream_sample_set.py`](../downstream_sample_set.py) (runner:
[`../scripts/build_downstream_sample_set.py`](../scripts/build_downstream_sample_set.py))
combines two exclusion sources into the final tables: the Stage 06 QC decisions
(`ignored_downstream = yes`) and the Stage 00 missing-mate exclusions. Included
samples get `downstream_cpDNA_use = include` and `downstream_mtDNA_use = include`;
excluded samples get an `exclusion_stage` (`step5_downstream_qc` or
`step0_manifest`), a reason, and evidence.

The stage guards its own headline number:

```python
DEFAULT_EXPECTED_INCLUDED = 275
...
if expected_included is not None and len(included_rows) != expected_included:
    raise DownstreamSampleSetError(
        f"Expected {expected_included} downstream included samples, found {len(included_rows)}.")
```

So if the arithmetic ever drifts — a sample added, a QC decision changed — the
stage *fails loudly* rather than silently producing a 274- or 276-sample set. A
test builds three analysis samples (one QC-dropped) plus one missing-mate sample
and asserts the two survivors are included while the QC-dropped and missing-mate
samples are excluded with the correct stages and reasons. `[TEST]`

This is where 278 complete pairs become 275: minus 2 missing mates (already
excluded at Stage 00, re-recorded here) and minus 3 low-coverage samples (Stage
06). `analysis_samples.tsv` had already dropped the missing mates, so the
275 = 278 − 3 arithmetic at this stage counts only the QC exclusions among the
included candidates. Every later stage reads `included_samples.tsv`.

## 9.5 The Python concepts here

- **Complement/merge interval algebra** (`complement_regions`,
  `merge_coordinate_pairs`) built from sorted coordinate pairs.
- **Coordinate conversion** in both directions ([Chapter 5,
  §5.5](./05-bioinformatics-file-formats.md)).
- **Cross-module import** of an entire stage's functions (06 from 02).
- **A hard invariant check** (`expected_included`) that turns a silent miscount
  into a raised error — a small idea with big consequences for trust.
- **Dict-of-counters** accumulation in `parse_track_depth_file`.

## 9.6 Failure modes

- **Malformed BED** (negative start, end ≤ start, too few fields) →
  `read_track_regions` raises `AlignmentError`. `[CODE]`
- **Track interval past the reference end** → `validate_track_regions` raises.
  `[CODE]`
- **Wrong number of IR intervals** → `MaskDefinitionError`. `[CODE]`
- **QC decision names a sample absent from the analysis set** →
  `DownstreamSampleSetError`, catching a stale or typo'd sample ID. `[CODE]`
- **Included count ≠ 275** → `DownstreamSampleSetError`. `[CODE]`
- **Biological failure: using the wrong track.** Nothing crashes if you
  hand-edit a config to call variants on the full mtDNA reference — but you would
  reintroduce the repeat artifacts Stage 03 removed. The track contract exists to
  make that mistake hard, not impossible. `[BIO]`

## 9.7 Exercises

1. **Trace.** A 30-bp reference has masked interval `10–20` (1-based inclusive).
   What population-site regions does `complement_regions` produce? Give them in
   1-based inclusive coordinates.
2. **Predict.** `interval_to_bed_fields("mitochondria", 4, 8, "x")` returns what?
   How many bases does the resulting BED interval cover?
3. **Predict.** You lower `DEFAULT_MTDNA_HIGH_CONFIDENCE_THRESHOLD` from 12 to
   10. Directionally, does the mtDNA population track get larger or smaller, and
   what is the trade-off for downstream mtDNA variants?
4. **Modify.** You want a third cpDNA track that is the small single-copy region
   only. Sketch the `Region` list and the `analysis_tracks.tsv` row you would
   add, and which downstream readers would need to opt into it.
5. **Debug.** Stage 07 raises `Expected 275 downstream included samples, found
   276`. Name two files you would diff to find the extra sample and the two
   exclusion sources that feed the count.
6. **Interpret.** Stage 06 says a sample has `cpdna_full_coverage` breadth 0.97
   but `cpdna_population_sites` breadth 0.99. Why can the population track's
   breadth be *higher* than the full-coverage track's?

Solutions in [Chapter 19](./19-solutions.md).

> Next: [Chapter 10 — From Reads to Alignments (Stages 08–11)](./10-variants-to-alignments.md)
