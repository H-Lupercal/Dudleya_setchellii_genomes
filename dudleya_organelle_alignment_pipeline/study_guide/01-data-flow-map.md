# Chapter 1 — The Data-Flow Map

> Part 1 of 4 · Foundations · Next: [Python
> Essentials](./02-python-essentials.md)

The fastest way to stop feeling lost in this repository is to see it as one
chain: each stage reads specific files, writes specific files, and hands the
next stage exactly what it needs. Nothing is magic and nothing is global. If you
can name the file that comes out of a stage, you can find the stage that reads
it next.

## 1.1 The whole pipeline as a provenance chain

```text
FASTQ files + population metadata
        |
        v
manifest and reference preflight (00-01)
        |
        v
pilot mapping and reference investigations (02-04)
        |
        v
masks, all-sample alignment, sample QC (05-07)
        |
        v
haploid calls, filters, SNP/callable alignments (08-11)
        |
        +----------------------+----------------------+------------------+
        v                      v                      v                  v
phylogenetic trees       PCA and clustering      population Fst     tool audit
(12, 14, 19, 20)         (15, 16, 18)            (17)               (13)
```

Read that top to bottom. The left column is a strict pipeline: you cannot call
variants (08) before you have masks (05) and alignments (06). The bottom row is
four independent consumers of the same upstream products. The tool audit (13)
sits off to the side because it depends on nothing but your installed software.

Two rules hold at *every* arrow, and they are the spine of the whole design:

1. **cpDNA and mtDNA are processed separately.** They share one reference FASTA
   and one BAM per sample, but every downstream product — variants, alignments,
   trees, PCA, clustering, Fst — is computed once for cpDNA and once for mtDNA.
   The biological reason is in [Chapter 6](./06-organelle-biology.md).
2. **QC tracks and population-genetic tracks are not interchangeable.** A region
   good enough to *measure a sample's coverage* is not necessarily a region you
   trust to *call a variant*. Stage 05 makes that distinction concrete as BED
   files. Stage 08 selects the intended population tracks by fixed track IDs, but
   the manifest's `purpose` field is not itself enforced; preserving the
   distinction still depends on correct track definitions ([Chapter
   9](./09-masks-alignment-and-sample-qc.md)).

## 1.2 The stage table

Each stage is identified by its `results/NN_.../` directory number — that number
is the canonical name used in every path. The authoritative version of this
table is [`../PROCESS.md`](../PROCESS.md); this copy adds the module and the
book chapter.

| Dir | Purpose | Runner script | Module | Chapter |
|---|---|---|---|---|
| 00 | Sample manifest + preflight | `scripts/build_sample_manifest.py` | [`manifest.py`](../manifest.py) | [7](./07-manifest-and-reference-preflight.md) |
| 01 | Reference validation, indexing, pilot selection | `scripts/prepare_reference_and_pilot.py` | [`prepare_reference_and_pilot.py`](../prepare_reference_and_pilot.py) | [7](./07-manifest-and-reference-preflight.md) |
| 02 | Pilot alignment (15 samples) | `scripts/run_pilot_alignment.py` | [`pilot_alignment.py`](../pilot_alignment.py) | [8](./08-pilot-mapping-and-investigations.md) |
| 03 | mtDNA repeat/placement investigation | investigation evidence | (consumed by 05) | [8](./08-pilot-mapping-and-investigations.md) |
| 04 | cpDNA inverted-repeat verification | investigation evidence | (consumed by 05) | [8](./08-pilot-mapping-and-investigations.md) |
| 05 | Analysis masks and tracks | `scripts/build_analysis_masks.py` | [`analysis_masks.py`](../analysis_masks.py) | [9](./09-masks-alignment-and-sample-qc.md) |
| 06 | All-sample alignment + track QC | `scripts/run_all_sample_alignment.py` | [`all_sample_alignment.py`](../all_sample_alignment.py) | [9](./09-masks-alignment-and-sample-qc.md) |
| 07 | Downstream include/exclude set (275) | `scripts/build_downstream_sample_set.py` | [`downstream_sample_set.py`](../downstream_sample_set.py) | [9](./09-masks-alignment-and-sample-qc.md) |
| 08 | Haploid variant calling | `scripts/run_variant_calling.py` | [`variant_calling.py`](../variant_calling.py) | [10](./10-variants-to-alignments.md) |
| 09 | Variant filtering / callable sites | `scripts/run_variant_filtering.py` | [`variant_filtering.py`](../variant_filtering.py) | [10](./10-variants-to-alignments.md) |
| 10 | Filtered SNP alignments | `scripts/run_snp_alignment.py` | [`snp_alignment.py`](../snp_alignment.py) | [10](./10-variants-to-alignments.md) |
| 11 | Full callable-site consensus alignments | `scripts/run_callable_consensus.py` | [`callable_consensus.py`](../callable_consensus.py) | [10](./10-variants-to-alignments.md) |
| 12 | ML trees, initial (no bootstrap) | `scripts/run_phylogenetic_tree.py` | [`phylogenetic_tree.py`](../phylogenetic_tree.py) | [11](./11-phylogenetic-trees.md) |
| 13 | Bioinformatics tool audit | `scripts/run_tool_audit.py` | [`tool_audit.py`](../tool_audit.py) | [14](./14-tool-audit.md) |
| 14 | Tree figures for Stage 12 | `scripts/run_tree_visualization.py` | [`tree_visualization.py`](../tree_visualization.py) | [11](./11-phylogenetic-trees.md) |
| 15 | PCA (cpDNA + mtDNA) | `scripts/run_pca_analysis.py` | [`pca_analysis.py`](../pca_analysis.py) | [12](./12-pca-and-clustering.md) |
| 16 | ADMIXTURE clustering, single-run | `scripts/run_admixture_analysis.py` | [`admixture_analysis.py`](../admixture_analysis.py) | [12](./12-pca-and-clustering.md) |
| 17 | Pairwise Fst + population summaries | `scripts/run_population_genetics.py` | [`population_genetics.py`](../population_genetics.py) | [13](./13-population-fst.md) |
| 18 | ADMIXTURE, five-replicate (final) | `scripts/run_admixture_analysis.py --replicates 5` | [`admixture_analysis.py`](../admixture_analysis.py) | [12](./12-pca-and-clustering.md) |
| 19 | ML trees, 1,000 UFBoot (final) | `scripts/run_phylogenetic_tree.py --bootstrap-replicates 1000` | [`phylogenetic_tree.py`](../phylogenetic_tree.py) | [11](./11-phylogenetic-trees.md) |
| 20 | Tree figures for Stage 19 (final) | `scripts/run_tree_visualization.py` | [`tree_visualization.py`](../tree_visualization.py) | [11](./11-phylogenetic-trees.md) |

Notice that stages 12/19, 14/20, and 16/18 reuse the *same module* with
different arguments. There is no separate "bootstrap tree" code — the same
`phylogenetic_tree.py` runs a fast search for Stage 12 and a 1,000-replicate
search for Stage 19. That is the initial-vs-final distinction, and it is worth
its own section.

## 1.3 Initial versus final: three analyses were run twice

Three analyses have a quick first pass and a rigorous final pass. The final pass
supersedes the first for interpretation and for the deliverables. When you read
a result, always know which one you are looking at. `[CODE]`

| First pass | Final pass | What changed |
|---|---|---|
| Stage 12 (ML trees) | Stage 19 (ML trees) | Final adds 1,000 ultrafast bootstrap replicates with BNNI. |
| Stage 14 (tree figures) | Stage 20 (tree figures) | Final renders the bootstrap-supported Stage 19 trees. |
| Stage 16 (ADMIXTURE) | Stage 18 (ADMIXTURE) | Final uses five seeded replicates per K and selects K by lowest *mean* CV error. |

Why keep both? Because the fast pass is cheap and lets you see the topology or
the clustering shape before committing hours of compute to bootstrapping. The
first passes are kept in the repo as historical record, not as deliverables. The
final deliverables and their exact paths are listed in
[`../results/final_deliverables_manifest.tsv`](../results/final_deliverables_manifest.tsv),
and the integrated narrative is
[`../results/organelle_population_report.md`](../results/organelle_population_report.md).

Any time this book cites a tree or a clustering result, it points you at the
*final* stage (19/20/18) unless it is explicitly discussing the first pass.

## 1.4 Run labels: `primary` and `smoke`

The completed analysis was run under the `primary` run label on the 275-sample
downstream set. That label is why you see `primary.` in so many filenames:
`primary.variant_calling_report.md`, `cpDNA.primary.callable_consensus.fa`, and
so on. `[CODE]`

There is also a 5-sample `smoke` run label. It exists only to validate the
variant-calling stage on a controlled handful of samples before spending compute
on all 275. It is **not** part of the primary results. When you see `smoke` in a
path, you are looking at a test run, not a deliverable. The mechanics of how a
run label threads through every output filename are in
[Chapter 3, §3.4](./03-reusable-code-patterns.md).

## 1.5 The numbers that flow through the chain

It helps to carry a few real numbers in your head as you read. These are
observed in the checked-in `primary` outputs. `[RESULT]`

| Checkpoint | cpDNA | mtDNA | Source |
|---|---:|---:|---|
| Samples discovered (all) | 280 | 280 | integrated report |
| Complete paired-end | 278 | 278 | integrated report |
| Downstream analysis set | 275 | 275 | Stage 07 |
| Callable-site alignment length | 124,538 | 44,930 | Stage 11 |
| Raw variant records | 2,475 | 190 | Stage 08 |
| Filtered SNPs | 2,015 | 146 | Stage 09 |
| Populations with codes (Fst) | 34 | 34 | Stage 17 |
| Lowest tested mean-CV K (five-replicate) | 8 | 8 | Stage 18 |

The story those numbers tell: the **analyzed cpDNA track** is much larger than the
trusted mtDNA track and contains roughly 14× more usable SNPs (2,015 vs 146).
(The full mitochondrial reference is actually larger: 243,359 bp versus 150,274
bp.) The smaller mtDNA marker set gives its tree, PCA, clustering, and Fst less
site-level resolution than their cpDNA counterparts; it does not by itself make
the cpDNA history biologically truer ([Chapter 6](./06-organelle-biology.md)).

## 1.6 How to navigate from here

- If Python syntax is what slows you down, read [Chapter
  2](./02-python-essentials.md) and [Chapter 3](./03-reusable-code-patterns.md)
  next; every pipeline chapter assumes them.
- If you want to understand a specific stage right now, jump to its walkthrough
  chapter from the table in §1.2.
- If you want to audit coverage — "is every module and function actually
  explained?" — use the [Module and Function
  Index](./21-module-and-function-index.md) in Part 4.

> Next: [Chapter 2 — Python Essentials for This Codebase](./02-python-essentials.md)
