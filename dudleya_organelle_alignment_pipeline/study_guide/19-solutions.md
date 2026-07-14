# Chapter 19 — Solutions

> Part 4 of 4 · Practice and Reference · Prev:
> [Exercises](./18-exercises.md) · Next: [Glossary](./20-glossary.md)

Solutions are keyed by chapter. Each explains the reasoning and points back to the
source or test. Where a question asks you to *modify* code, the answer describes
the change; it does not ask you to run anything.

## Chapter 7 — Manifest and Reference Preflight

**7.1 Trace.** `parse_fastq_path` splits off the `_S77_L006_R2_001.fastq.gz`
suffix, leaving prefix `ABMU_HOR_LP_140_Du-410`. `MAIN_STANDARD_RE` matches:
popcode `ABMU_HOR` (greedy up to `_LP_`), `lp_id=LP_140`, `du_id=Du-410`, so
`naming_profile=main_standard`, `read=R2`. If the popcode is absent from the CSV,
`infer_species_from_popcode("ABMU_HOR")` — not `CY_`, but starts with `AB` — gives
**D. abramsii**, with `metadata_status = popcode_not_in_csv`.

**7.2 Predict.** Two R1 and two R2: `determine_pair_status` sees equal, >1 counts →
`complete_multi_file`; `determine_analysis_status` maps that to
`review_before_primary_analysis`. It is therefore not in `analysis_samples.tsv`,
so it never reaches the pilot; even if forced in, `read_alignment_samples` requires
`pair_status == "complete"` and would skip it.

**7.3 Predict.** With the CSV deleted, `load_population_codes` returns `{}` (it
checks `csv_path.exists()`). For `CY_RED_...`, `population` is `None`, so species
falls back to `infer_species_from_popcode("CY_RED")` = **D. cymosa**, and
`metadata_status = popcode_not_in_csv`. No crash.

**7.4 Modify.** Add e.g. `INITIAL_DUSE_RE = re.compile(r"^(?P<du_id>DUSE_\d+)$")`
and a branch in `classify_prefix` returning `("initial_duse", "", du_id, "")`,
placed **before** the final `unrecognized` return. Order is safe because
`MAIN_STANDARD_RE` needs `_LP_\d+_Du-` and will not match `DUSE_412`.

**7.5 Debug.** (a) Filename parse — check `pairing_report.tsv` for an
`unparsed_fastq_name` issue. (b) Batch grouping — `infer_batch` keys samples by
batch, so the same ID in two batches splits into two rows; check `samples.tsv`.
(c) Pair status — check the sample's `pair_status`/`analysis_status` in
`samples.tsv`; only `complete` reaches `analysis_samples.tsv`.

**7.6 Interpret.** Not a Stage 00 bug. Stage 00 has aligned nothing, so it cannot
know coverage. The three missing samples were dropped by **Stage 06** coverage QC
and formalized in **Stage 07** (278 complete − 3 low-coverage = 275). Pairing is
necessary, not sufficient.

## Chapter 8 — Pilot Mapping and Investigations

**8.1 Trace.** `parse_idxstats_file` skips the `*` row →
`{"chloroplast": 5000, "mitochondria": 40}`. In `build_sample_summary`, mt = 40
< 100 fires `low_mitochondria_mapped_reads` (plus any low-breadth note depending
on the depth file). cp = 5000 is fine.

**8.2 Predict.** BAM exists and `--force` is off, so the alignment is reused;
`--refresh-qc` regenerates index, flagstat, idxstats, depth. `commands.tsv`
records steps `index`, `flagstat`, `idxstats`, `depth`. The BAM is not rebuilt.

**8.3 Predict.** Depths 4 and 6 over a 2-base reference: `total_depth = 10`,
`mean_depth = 5.0`, `breadth_ge_1x = 2/2 = 1.0`, `breadth_ge_5x = 1/2 = 0.5`
(only the 6 clears 5).

**8.4 Modify.** Add `bases_ge_20x` to `OrganelleMetrics` plus a `breadth_ge_20x`
`@property`; add the counter and `if depth >= 20` line in `parse_depth_file`; add
the two fields in `build_organelle_summary_rows` (and, if you want it in the
per-sample row, `build_sample_summary`).

**8.5 Debug.** cp 0.999 vs mt 0.02: this is a low-mtDNA-coverage sample. It is
plausibly biological (little mtDNA in the library) rather than a mapping failure,
because cpDNA mapped fine against the same reference. Confirm by checking the
mtDNA mapped-read count and the mtDNA investigation before deciding; the QC note
`low_mitochondria_breadth_ge_1x` will already be set.

**8.6 Interpret.** The 0.96 permissive median counts reads at MAPQ ≥0, which
repeats inflate. The 44,930 bp is the high-MAPQ, uniquely-placed subset supported
in ≥12 pilot samples. Demanding unique placement collapses much of the apparent
coverage — hence the small trusted track.

## Chapter 9 — Masks, Alignment, and Sample QC

**9.1 Trace.** Mask `10–20` on a 30-bp reference: `complement_regions` emits
`1–9` (before the mask) and `21–30` (after), in 1-based inclusive coordinates.

**9.2 Predict.** `interval_to_bed_fields("mitochondria", 4, 8, "x")` →
`["mitochondria", "3", "8", "x"]`; the BED interval `3–8` (half-open) covers
`8 − 3 = 5` bases.

**9.3 Predict.** Lowering the threshold 12 → 10 selects the intervals recorded at
the 10-sample support level, which are generally more/larger, so the mtDNA
population track **grows**. Trade-off: more mtDNA sites but lower confidence
(regions unique in fewer samples), risking the repeat artifacts the track exists
to avoid.

**9.4 Modify.** Build a `Region` list for the small single-copy region, add a row
to `analysis_tracks.tsv` (e.g. `cpdna_small_single_copy`, purpose `annotation`),
and opt any consumer in by adding its `track_id` to that consumer's wanted set —
`read_variant_tracks` in `variant_calling.py` and `read_consensus_inputs` in
`callable_consensus.py` both hard-code the two population track IDs. Nothing uses
a new track until a reader names it.

**9.5 Debug.** Found 276 instead of 275: diff `included_samples.tsv` against the
expected set, and check the two exclusion sources — the QC decisions
(`ignored_downstream = yes` in `downstream_sample_qc_decisions.tsv`) and the
missing-mate rows in `excluded_samples.tsv`. An extra included sample means one
expected QC exclusion is missing.

**9.6 Interpret.** The population track omits the duplicate IR copy. If that
region is lower-covered in a sample, removing it from the denominator raises the
computed breadth, so the population-track breadth can exceed the full-reference
breadth.

## Chapter 10 — From Reads to Alignments

**10.1 Trace.** `T C` with genotypes `1 0 .` → sample 1 = ALT = **C**, sample 2 =
REF = **T**, sample 3 = **N**.

**10.2 Predict.** Raw sites {3, 6}, filter keeps {3} → failed = {6}. Every sample
gets **N** at position 6, because the caller saw a variant there but it failed
filtering, so neither the reference nor a called allele is trusted.

**10.3 Predict.** Raising `--min-depth` 1 → 5 masks more positions → `missing_bases`
**increases** and the number of informative (non-`N`) columns **decreases**.

**10.4 Modify.** Change `DEFAULT_MAX_MISSING_FRACTION` to `0.3` and the
`-i 'F_MISSING<=...'` flag. Allowing more missingness lets more sites pass, so the
filtered SNP count goes **up**.

**10.5 Debug.** "Raw VCF sample order does not match included samples": the VCF was
built from a different or older sample set, or `included_samples.tsv` changed
after calling. Compare the VCF `#CHROM` sample columns to the `sample_id` order in
`included_samples.tsv`.

**10.6 Interpret.** cpDNA ≈ 2,015 / 124,538 ≈ 0.0162 SNP/bp (1 per ~62 bp); mtDNA
≈ 146 / 44,930 ≈ 0.0033 (1 per ~308 bp). cpDNA is ~5× denser, so the mtDNA tree
has far less signal per tip and lower expected support.

## Chapter 11 — Phylogenetic Trees

**11.1 Trace.** `fast = not False and not 1000` = `True and False` = **False**;
`build_iqtree_command` appends `-B 1000 --bnni` and **not** `--fast`.

**11.2 Predict.** `--full-search` → `fast = not True and ...` = False, bootstrap
default 0 → `method = "iqtree_ml"`, and `--fast` is absent.

**11.3 Predict.** `compute_tree_figure_size(5)` → height `max(6.0, min(80.0,
2.8)) = 6.0` → `(14.0, 6.0)`. `(275)` → `275*0.16+2 = 46.0` → `(14.0, 46.0)`. The
`min(80.0, ...)` cap keeps a very large tree from producing an unusably tall
figure.

**11.4 Modify.** `run_phylogenetic_tree.py --bootstrap-replicates 5000
--output-dir .../21_...` then point `run_tree_visualization.py --tree-dir` at it.
No module changes — `--bootstrap-replicates` is already a parameter.

**11.5 Debug.** The post-run guard `if not treefile_path.exists() or size == 0`
raises even on exit 0. Inspect `<prefix>.log`/`<prefix>.iqtree`; the treefile is
missing or empty (wrong prefix/output-dir, or an alignment IQ-TREE could not use).

**11.6 Interpret.** Before calling it biological discordance, examine mtDNA
resolution (146 SNPs and UFBoot on the conflicting branch) plus sample-specific
callability. Overall mtDNA missingness is only 0.2534%, so do not describe the
matrix as heavily missing; use Stage 11 to check whether missingness is
concentrated in the samples involved.

## Chapter 12 — PCA and Clustering

**12.1 Trace.** Column `A, A, T, N`: alleles `{A, T}` → `A=0.0`, `T=1.0`. Encoded
`0, 0, 1, nan`; `nanmean = 1/3 ≈ 0.333` replaces the `N`.

**12.2 Predict.** `G,G,G,G` is monomorphic (`len(alleles) < 2`) → dropped.
`A,A,N,N` has only one observed allele → also dropped.

**12.3 Predict.** `["admixture", "--cv", "--seed=100", "-j8", "mtDNA.ped", "5"]`
(the command uses `genotype_path.name`).

**12.4 Modify.** Change `choose_plot_group` in `pca_analysis.py`. With species and
popcode both blank it already returns `naming_profile or "unresolved"`, so
returning `naming_profile or "unresolved"` unconditionally gives profile-only
coloring.

**12.5 Debug.** If mean CV error is lowest at K=1 and rises with K, K=1 is the
best prediction among the tested values under this model. It does not prove an
absence of biological structure, especially with pseudo-diploid linked markers.

**12.6 Interpret.** A monotonic CV decline to the maximum tested K means the sweep
did not bracket an interior optimum. Under the additional pseudo-diploid and
linked-marker model violations, the honest reading is “K=8 had the lowest CV
error among K=1–8,” not “at least” or “exactly” eight biological groups.

## Chapter 13 — Population Fst

**13.1 Trace.** Site 1: A = `{A:2}` → H=0; B = `{T:2}` → H=0; combined `{A:2,T:2}`
→ H_T = 1 − (0.25+0.25) = 0.5; H_S = 0. Contribution 0.5/0.5. Same at site 2 →
Fst = 1.0/1.0 = **1.0**, informative_sites = 2, matching the test.

**13.2 Predict.** `AA, AT, TT` → 3 distinct haplotypes; diversity =
`(3/2)·(1 − 3·(1/3)²) = (3/2)·(2/3) = 1.0`.

**13.3 Predict.** 34 populations → `34·33/2 = 561`; 10 populations → `45`.

**13.4 Modify.** After `group_sequences_by_population`, filter
`groups = {p: r for p, r in groups.items() if len(r) >= 3}`. Fewer populations →
fewer `combinations`, so the pairwise comparison count drops.

**13.5 Debug.** mtDNA Fst 0.0 where you expect a difference: likely missingness
(sites skipped because a population has no data) or monomorphic combined sites.
The `informative_sites` column in the pairwise table tells you how many sites
actually contributed; a low count means "no resolution," not "no difference."

**13.6 Interpret.** cpDNA 0.4 vs mtDNA 0.0 for one pair: first check whether
mtDNA's 146 SNPs leave too few informative sites for that pair. Overall mtDNA
missingness is low, but concentrated missing data could still reduce the
`informative_sites` count. A remaining difference is discordance between two
cytoplasmic markers, not automatically independent or specifically maternal
evidence.

## Chapter 14 — Tool Audit

**14.1 Trace.** `plink` resolves to `None`, so the loop tries `plink2`, finds it,
runs the version check, and returns `status=FOUND`, `executable="plink2"`.

**14.2 Predict.** `python3` present but no scikit-learn: the resolver finds
`python3`, the runner raises on import, the `except` returns `status=MISSING`,
`path=/…/python3`, note `"Executable found, but version/import check failed: ..."`.

**14.3 Predict.** All `required_current` found, `admixture` missing:
`ready_for_current_pipeline = True`, `ready_for_remaining_goal = False`
(`admixture` is in `missing_required_remaining`).

**14.4 Modify.** `ToolSpec("mafft", ("mafft",), "recommended_remaining",
"multiple sequence alignment", ("--version",))` added to the `TOOL_SPECS` tuple.
`recommended_remaining` because no current stage uses it.

**14.5 Debug.** The spec already lists `("iqtree", "iqtree2")`, so a MISSING
report means neither is on `PATH` in the audit's environment — a `PATH` problem.
Run with the `.tools/bioconda-env/bin` prefix; do not change the spec or resolver.

**14.6 Interpret.** A bare `python3` says nothing about installed libraries, but
Stage 15 (PCA) imports `sklearn` and `matplotlib`. Trusting `python3` alone would
pass the audit yet the PCA stage would crash at import — so the audit imports the
libraries to check them.

## Chapter 18 — Warm-up, Integrative, and Test exercises

**W1.** 16 test files, 70 tests, running in about a second. The heavy binaries
(`bwa`/`iqtree`/`admixture`/`plink`) are never invoked: command-builder tests
assert argument *lists* without executing tools, and pure-logic tests use tiny
fixtures. Non-standard Python dependencies are still required by rendering
tests: the PCA test uses NumPy, scikit-learn, and matplotlib, while
`test_render_tree_figure_writes_png_pdf_and_svg` and
`test_write_tree_visualization_outputs_records_summary_and_report` use Biopython
and matplotlib. Run the suite in the pipeline environment—the “import to check”
point of [Chapter 14](./14-tool-audit.md).

**W2.** `test_build_callable_consensus_applies_depth_variants_and_failed_site_mask`
asserts `"CTTNNG"`/`"CGTANG"`; it pins Stage 11 (callable consensus).

**I1.** 2,475 raw (Stage 08) → 2,015 filtered (Stage 09) → 2,015 SNP columns
(Stage 10). The 460 removed are multiallelic, indel, singleton, or too-missing.
The checked-in summaries do not attribute the drop per flag, so the honest way to
find the biggest contributor is to rerun Stage 09 changing one flag at a time —
in organelle SNP sets, the singleton filter (`--min-ac 2:minor`) and the
missingness filter usually remove the most.

**I2.** More permissive filters (`--max-missing-fraction 0.5`,
`--min-minor-allele-count 1`) let more sites pass: cpDNA filtered SNP count **up**,
Stage 10 width **up**, Stage 15 `retained_sites` **up** (minus any newly
monomorphic columns PCA drops). No stage raises; the cost is more low-information,
noisier sites.

**I3.** Stage 08 with `--run-label smoke` writes `smoke.variant_calling_summary.tsv`.
Stage 09 defaulting to `primary` looks for `primary.variant_calling_summary.tsv`,
does not find it, and raises `VariantFilteringError`. The mismatch is explained by
`labeled_output_name`, which prefixes every output with the run label
([Chapter 3, §3.4](./03-reusable-code-patterns.md)).

**I4.** Calling on `cpdna_full_coverage` reintroduces the duplicate IR copy (and
lower-MAPQ regions), so IR variants are double-counted and the cpDNA raw/filtered
SNP counts rise with spurious, correlated sites. It violates the "QC tracks are
not population tracks" rule ([Chapters 1](./01-data-flow-map.md) and
[9](./09-masks-alignment-and-sample-qc.md)).

**I5.** In `pilot_alignment.py`, touch `OrganelleMetrics` (add `bases_ge_20x` and
a `breadth_ge_20x` property), `parse_depth_file` (count ≥20×), and
`build_organelle_summary_rows` (emit the fields). Stage 06 reuses alignment
execution from that module, but its track-aware depth summaries are separate. In
`all_sample_alignment.py`, also update `TrackMetrics`,
`initialize_track_counters`, `parse_track_depth_file`, and
`build_track_summary_rows`. The ≥20× metric does **not** flow into Stage 06
automatically ([Chapter 3, §3.7](./03-reusable-code-patterns.md)).

**I6.** The two breadths reconcile via the duplicate-IR exclusion (§9.6). But
"mostly `N`" is **inconsistent** with a population-track breadth of 0.99 at
`min_depth = 1`: if 99% of positions have ≥1× coverage, the consensus is mostly
reference bases, not `N`. Suspect a mismatched or wrong depth file, a raised
`--min-depth`, or a depth/consensus file from a different run. The checked-in
mtDNA matrix is also not mostly `N`: its overall missingness is 0.2534%.

**I7.** A defensible sentence: "For cpDNA, populations P and Q separate on PC1,
fall in distinct UFBoot-98 clades, and show high relative pairwise Fst (0.45)
`[RESULT]`, consistent with differentiated *cpDNA haplotypes*
`[BIO]`; the mtDNA marker does not resolve this split (overlapping PCA, low tree
support, Fst 0.02), and none of this establishes nuclear divergence or species
status."

**I8.** S1 = `"CTTNNG"` — traced fully in [Chapter 23](./23-capstone-sample-trace.md),
matching `test_build_callable_consensus_...`.

**T1.** IR copies `5–8` and `15–18` on a 20-bp reference → `cpdna_population_sites`
= `1–14` and `19–20` (everything except the masked later copy). Confirm against
the test's assertion.

**T2.** `-m2 -M2` (biallelic — clean, interpretable sites); `-v snps` (SNPs only —
drop indels); `--min-ac 2:minor` (no singletons — likely errors); `F_MISSING<=0.2`
(≤20% missing — enough samples genotyped). See [Chapter 10,
§10.4](./10-variants-to-alignments.md).

**T3.** K=3 has the lower mean (0.21 vs 0.31), so `is_best_mean_k = yes` for K=3.
To make K=2 win, set its replicate CV errors below K=3's (e.g. K=2 → 0.20, 0.20;
K=3 → 0.30, 0.30); then the assertion flips to K=2.

**T4.** Open-ended: any small fixture + single call + single assertion that pins a
behavior not already covered (e.g. that `median_breadth` averages the middle two
values for an even count) is a valid answer, as long as you can state the
behavior it protects.

> Next: [Chapter 20 — Glossary](./20-glossary.md)
