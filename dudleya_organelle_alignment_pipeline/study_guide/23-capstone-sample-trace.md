# Chapter 23 — Capstone: Tracing One Sample End to End

> Part 4 of 4 · Practice and Reference · Prev:
> [External Tool and Format Reference](./22-external-tool-and-format-reference.md)

This capstone follows a single sample the whole length of the pipeline, then
hand-traces the one genuinely intricate computation — the callable-site consensus
— on the exact six-base example the test pins. If you can follow both, you
understand the pipeline.

## 23.1 One sample, conceptually, from FASTQ to interpretation

Take `CY_ALA_LP_298-Du-767`, a *D. cymosa* sample with the tricky
hyphen-separated name. Watch what each stage does *to this one sample*.

**Stage 00 — identity.** Its two FASTQs
(`CY_ALA_LP_298-Du-767_S289_L005_R1_001.fastq.gz` and the `R2`) are parsed:
`FASTQ_NAME_RE` peels the suffix, `MAIN_STANDARD_RE` reads
popcode `CY_ALA`, `lp_id LP_298`, `du_id Du-767`. Exactly one R1 and one R2 →
`pair_status = complete` → `analysis_status = include_primary_paired_end`. The
popcode resolves in the CSV (or infers *D. cymosa* from `CY_`). It lands in
`analysis_samples.tsv`. ([Chapter 7](./07-manifest-and-reference-preflight.md))

**Stage 01 — maybe a pilot pick.** If it is diverse enough to be selected, it
joins the 15-sample pilot; otherwise it waits for the all-sample run. Either way
the reference is validated and indexed first. ([Chapter 7](./07-manifest-and-reference-preflight.md))

**Stage 06 — mapping and QC.** `bwa mem | samtools view -F 4 | samtools sort`
produces `CY_ALA_LP_298-Du-767.organelle.sorted.bam` (note the safe name keeps
the hyphen, which is filename-legal). `idxstats` splits its reads into cpDNA and
mtDNA; `depth` measures coverage; the track summary reports its breadth on
`cpdna_population_sites` and `mtdna_high_confidence_unique`. Suppose its breadth
clears the screens — it is *not* one of the three excluded low-input samples.
([Chapter 9](./09-masks-alignment-and-sample-qc.md))

**Stage 07 — it makes the cut.** With `downstream_cpDNA_use = include` and
`downstream_mtDNA_use = include`, it is one of the 275. ([Chapter 9](./09-masks-alignment-and-sample-qc.md))

**Stages 08–11 — it becomes columns.** Its BAM joins the 275-sample `mpileup`.
At each variant site it contributes one haploid allele (`--ploidy 1`). After
filtering, it is one row in the cpDNA SNP alignment (2,015 columns) and one row
in the mtDNA SNP alignment (146 columns), and one record in each callable-site
consensus (124,538 and 44,930 columns), with `N` wherever its own depth was too
low or a raw variant failed filtering. ([Chapter 10](./10-variants-to-alignments.md))

**Stages 12–20 — it becomes a point, a bar, a tip.** In the cpDNA PCA it is one
dot, colored `D. cymosa_CY_ALA`; in the ADMIXTURE plot one stacked bar sorted
next to its `CY_ALA` neighbors; in the Stage 19 tree one tip, `CY_ALA_LP_298-Du-767`,
with UFBoot support on the branches around it. In Fst it contributes to the
`CY_ALA` population's diversity and its pairwise comparisons — *if* its popcode
resolved. ([Chapters 11](./11-phylogenetic-trees.md)–[13](./13-population-fst.md))

**Interpretation — what its position means.** If it clusters tightly with other
`CY_ALA` samples across cpDNA PCA, tree, and Fst, that is consistent evidence that
their cpDNA sequences share recent cytoplasmic history. It is **not**
evidence about their nuclear ancestry, species status, or admixture — one linked
cytoplasmic locus cannot carry those claims ([Chapter 17](./17-uncertainty-bias-and-limits.md)).

## 23.2 The one computation worth doing by hand

The callable-site consensus ([Chapter 10, §10.6](./10-variants-to-alignments.md))
is the pipeline's most intricate per-sample logic. The test
`test_build_callable_consensus_applies_depth_variants_and_failed_site_mask` traces
it on a tiny example; here is that example, step by step, so the algorithm is
concrete.

### The inputs

- **Reference** `chloroplast`: `ACGTACGT` — positions 1–8 are
  `A C G T A C G T`.
- **BED track**: `chloroplast 1 7` (0-based half-open) → 1-based positions **2–7**
  → six columns.
- **Samples**: S1 and S2 (both included).
- **Filtered VCF** (survived filtering): position 3, `G→T`, S1 = `1` (ALT), S2 =
  `0` (REF).
- **Raw VCF** (all calls): position 3 (same as above) **and** position 6, `C→A`,
  `LowQual`, both samples `1`.
- **Depth files** (min_depth = 1): S1 covers positions 2, 3, 4, 6, 7 (**not 5**);
  S2 covers 2, 3, 4, 5, 6, 7 (all).

### Step 1 — the template

Reference bases at positions 2–7 are `C G T A C G`. Column indices 0–5 map to
positions 2–7:

| index | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| position | 2 | 3 | 4 | 5 | 6 | 7 |
| reference | C | G | T | A | C | G |

Both samples begin as `CGTACG`.

### Step 2 — the failed-site set

Raw variant positions inside the track: {3, 6}. Filtered positions: {3}. The
failed set is `{3, 6} − {3} = {6}`, which is **index 4**. Every sample will get
`N` at index 4, because a raw call there failed filtering — we trust neither the
variant nor the reference.

### Step 3 — overlay the filtered SNP

Filtered SNP at position 3 (index 1), `G→T`: S1 = `1` → ALT = **T**; S2 = `0` →
REF = **G** (unchanged). After this step:

- S1: `C T T A C G`
- S2: `C G T A C G`

### Step 4 — the masks (per sample)

Apply the failed-site mask (index 4 → `N`) and the depth mask (positions not
covered → `N`).

**S1**: index 4 → `N` (failed). Depth: position 5 = index 3 is **not** covered →
`N`. Result: `C T T N N G` = **`CTTNNG`**.

**S2**: index 4 → `N` (failed). Depth: all positions covered → no extra `N`.
Result: `C G T A N G` = **`CGTANG`**.

### Step 5 — the tallies

- `filtered_variant_sites` = 1 (position 3 applied).
- `masked_failed_variant_sites` = 1 (position 6).
- `missing_bases` = S1's two `N`s (indices 3, 4) + S2's one `N` (index 4) = **3**.

These are exactly the values the test asserts: `sequences["S1"] == "CTTNNG"`,
`sequences["S2"] == "CGTANG"`, `filtered_variant_sites == 1`,
`masked_failed_variant_sites == 1`, `missing_bases == 3`. `[TEST]`

### Why each `N` is there — the whole point

The two `N`s in S1 come from *two different reasons*, and telling them apart is
the pipeline's core discipline:

- Index 4 (position 6) is `N` for **both** samples because a called variant failed
  quality filtering — a *site-level* trust problem.
- Index 3 (position 5) is `N` for **S1 only** because S1 had no read coverage
  there — a *sample-level* callability problem.

The consensus records absence of evidence as `N`, and it distinguishes "we
couldn't call a good variant here" from "this sample wasn't covered here." Every
downstream method then treats those `N`s as missing data, not as a base — which
is why the interpretation chapters keep insisting that the width of an alignment
is not the amount of evidence in it.

## 23.3 You have finished the book

You can now, for any sample: parse its identity, follow its reads to a BAM, see
how QC keeps or drops it, watch it become haploid variant calls, SNP columns, and
a masked consensus record, and read its place in a tree, PCA, cluster, and Fst —
with the correct caution about what one largely linked cytoplasmic locus can and cannot tell
you.

Use the [Module and Function Index](./21-module-and-function-index.md) to jump to
any symbol, and reread [Chapter 17](./17-uncertainty-bias-and-limits.md) whenever
you are about to write down a conclusion. The pipeline's whole design — separate tracks, explicit masks,
loud failures, labeled evidence — exists so that the confident claims are actually
worth trusting. Read it in that spirit.
