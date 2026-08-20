# Chapter 17 — Uncertainty, Bias, and the Limits of These Claims

> Part 3 of 4 · Interpretation · Prev: [Reading PCA,
> Clustering, and Fst](./16-pca-clustering-fst-interpretation.md) · Next:
> [Exercises](./18-exercises.md)

This chapter collects every caution scattered through the book into one place, so
you can audit any claim against it. It is organized as a hierarchy: the one
limit that dominates everything, then the biological biases, then the
computational ones, then a claim-by-claim ledger.

## 17.1 The dominant limit: each organelle is one cytoplasmic locus

Everything this pipeline produces — trees, PCA, clustering, Fst, for both cpDNA
and mtDNA — describes **organelle-haplotype variation**. Each organelle is treated
as one largely linked cytoplasmic locus. Maternal inheritance is common in plants,
but this repository does not establish inheritance direction in *Dudleya*, and
plant mitochondrial DNA can recombine among repeats. No number of SNPs,
bootstrap replicates, or fitted components turns either organelle into a nuclear
genome scan. `[BIO]`

The consequences, which no downstream statistic can undo:

- **These analyses cannot establish nuclear admixture.** Cytoplasmic loci carry
  no direct measurement of genome-wide biparental nuclear ancestry. The
  pseudo-diploid ADMIXTURE analysis is an exploratory projection under violated
  ploidy and linkage assumptions, not a nuclear ancestry analysis.
- **A gene tree is not a species tree.** A well-supported organelle clade shows
  shared organelle history, which can differ from the species tree through
  lineage sorting, organelle capture, or cytoplasmic introgression.
- **"Sample size" is misleading.** 2,015 cpDNA SNPs are 2,015 points on *one*
  genealogy, not 2,015 independent loci. Confidence in the topology is not the
  same as confidence about the organism's history.
- **cpDNA and mtDNA are separate but not automatically independent markers.**
  Comparing them is informative, but shared transmission can couple their
  histories, and neither measures the nuclear genome.

If you remember one thing from this book, remember this section. Every other
caution is a detail beneath it.

## 17.2 Biological biases baked into the data

**Reference bias.** All samples map to one Dudleya organelle reference; reads too
divergent from it map poorly or not at all, so variation in highly divergent
regions is under-seen ([Chapter 6, §6.6](./06-organelle-biology.md)). This cannot
be removed with a single reference; the pipeline mitigates it by restricting
analysis to trustworthy regions and reporting breadth, not by pretending it is
absent. `[BIO]`

**mtDNA repeats and paralogy.** Plant mitochondria recombine among repeats and
carry sequence of plastid/nuclear origin. Reads from those regions map
ambiguously, so the pipeline keeps only ~44,930 bp of uniquely mappable mtDNA
([Chapter 8](./08-pilot-mapping-and-investigations.md)). The price is only 146
usable mtDNA SNPs — every mtDNA result is lower-resolution than its cpDNA
counterpart. `[RESULT]`

**Chloroplast inverted repeat.** Kept as a single copy to avoid double-counting
correlated variation ([Chapter 6, §6.3](./06-organelle-biology.md)); the
population alignment is 124,538 bp, not the full 150,274. `[RESULT]`

**Heteroplasmy and haploidy assumption.** The pipeline treats each organelle as
one haplotype (`--ploidy 1`). Real heteroplasmy or nuclear-integrated organelle
sequences (numts/nupts) would violate this, and such sites are handled as
missing/filtered rather than modeled. `[BIO]`

## 17.3 Computational and statistical limits

**Callability and missingness.** A base is `N` unless it had enough confident
coverage. Across the full matrices, mtDNA has 31,313 missing cells (**0.2534%**)
and cpDNA has 127,485 (**0.3722%**). `[RESULT]` Every method handles `N`
differently—PCA mean-imputes, Fst excludes missing observations, and ADMIXTURE
codes them `0 0`. Overall missingness is low, but concentrated sample- or
site-level missingness can still matter; inspect its distribution.

**The MAPQ asymmetry.** Callability is defined from permissive-MAPQ Stage 06
depth, while variants are called at MAPQ ≥20 ([Chapter 10,
§10.6](./10-variants-to-alignments.md)). "Callable" is therefore a slightly more
generous notion than "confidently variant-callable." `[CODE]`

**The K boundary and model mismatch.** ADMIXTURE's lowest tested mean CV error is
K=8, the top of the 1–8 range, with a monotonically decreasing curve. The sweep
did not bracket an interior optimum and supplies no lower bound on biological
group count. Pseudo-diploid encoding and linked organelle SNPs also violate core
model assumptions ([Chapter 16,
§16.2](./16-pca-clustering-fst-interpretation.md)). `[RESULT]`

**Filtering choices are thresholds, not truths.** Biallelic-only, SNP-only,
minor-allele-count ≥2, ≤20% missing (Stage 09) are defensible defaults, but they
are choices: they discard indels, singletons, and multiallelic sites. Different
thresholds would give different SNP counts and slightly different downstream
pictures. `[CODE]`

**Model and estimator choices.** Trees assume `GTR+F+G4`; Fst uses one Nei-style
estimator. Both are reasonable and both are choices; numbers are not portable
across models/estimators. `[CODE]`

**Software failure modes** (all raised loudly, never silent): missing tools,
malformed paths, failed subprocesses, partial outputs, missing mates, invalid
sample metadata, malformed TSV/FASTA/BED/BAM/VCF, sample-order mismatches. The
pipeline's design preference throughout is a loud crash over a plausible wrong
answer ([Chapter 2, §2.9](./02-python-essentials.md)).

## 17.4 The claim ledger

A compact audit of what the pipeline supports, by evidence type. Use it to grade
any statement you are tempted to write.

| Claim | Supported? | Strongest evidence |
|---|---|---|
| "The tests define the intended behavior of each function." | Yes | `[TEST]` — 16 test files |
| "cpDNA yields 2,015 filtered SNPs, mtDNA 146." | Yes | `[RESULT]` Stage 09 |
| "The pipeline calls variants haploid, per organelle, on trusted tracks." | Yes | `[CODE]` Stage 05/08 |
| "This cpDNA clade is well supported (UFBoot ≥95)." | Yes, as *organelle* history | `[RESULT]` + `[BIO]` |
| "These individuals form a species / are related overall." | No | one unrooted locus |
| "The ADMIXTURE plot shows biological ancestry/assignment proportions." | No | pseudo-diploid encoding + linked sites |
| "There are at least/exactly 8 biological groups." | No | boundary K + model mismatch |
| "mtDNA shows no structure, so there is none." | No | only 146 SNPs / limited resolution |
| "cpDNA Fst 0.4 equals a nuclear Fst of 0.4." | No | estimator/marker-specific |
| "These organelle analyses establish complete species history." | No | §17.1 |

## 17.5 How to write a defensible sentence

A good claim from this pipeline names three things: the **organelle**, the
**evidence type**, and the **limit**. For example:

> "For cpDNA, populations A and B fall in separate, UFBoot-≥95 clades in the
> Stage 19 tree and show high relative pairwise Fst `[RESULT]`, consistent with
> differentiated *cpDNA haplotypes* `[BIO]`; this does not by
> itself establish nuclear divergence or species status."

That sentence would survive review. A sentence that drops the organelle
qualifier, the evidence tag, or the limit would not. When in doubt, add the
qualifier — the biology in §17.1 makes it mandatory, not optional.

## 17.6 The point of all this caution

The pipeline is careful *on purpose*: separate tracks, explicit masks, loud
failures, labeled evidence, an integrated report that states its own caveats. The
caution is not timidity; it is what makes the confident claims — "here are 2,015
trustworthy cpDNA SNPs, here is a bootstrap-supported organelle tree" — actually
trustworthy. Reading the results well means matching your confidence to the
evidence type, every time.

> Part 3 complete. Next: [Chapter 18 — Exercises](./18-exercises.md)
