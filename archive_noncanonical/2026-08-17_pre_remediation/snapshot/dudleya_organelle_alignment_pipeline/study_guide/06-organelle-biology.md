# Chapter 6 — Organelle Biology That Changes the Code

> Part 1 of 4 · Foundations · Prev: [Bioinformatics
> File Formats](./05-bioinformatics-file-formats.md) · Next: [Manifest and
> Reference Preflight](./07-manifest-and-reference-preflight.md)

You know the biology of chloroplasts and mitochondria. This chapter is short on
definitions and long on *consequences*: the specific ways organelle genomes
force design decisions in this pipeline that a nuclear-genome pipeline would not
make. Every one of these ideas shows up later as a threshold, a mask, an
encoding, or a caveat.

## 6.1 Haploidy → `--ploidy 1` and single-allele genotypes

An organelle genome is present in many copies per cell, but they are (barring
heteroplasmy) effectively one haplotype: there is no pair of homologous
chromosomes, no heterozygotes, no phasing problem. Computationally, this is a
gift and a constraint.

- The variant caller runs with **`bcftools call --ploidy 1`**, so each sample
  gets a single allele per site, written as `0`, `1`, or `.` in the VCF, not
  `0/1` ([Chapter 4, §4.4](./04-shell-and-external-tools.md)).
- A genotype decodes to exactly one base, so an alignment column is a single
  letter per sample. There is no heterozygous `R`/`Y` ambiguity code anywhere in
  this pipeline.
- The diversity math is haplotype-based: `compute_haplotype_diversity` counts
  *unique whole sequences* as haplotypes ([Chapter 13](./13-population-fst.md)),
  which only makes sense because each sample *is* one haplotype.

The one place haploidy fights the tooling is ADMIXTURE, which assumes diploids.
The pipeline works around it by duplicating each haploid base into a homozygous
"pseudo-diploid" pair ([Chapter 5, §5.7](./05-bioinformatics-file-formats.md)).
That is a *tooling* accommodation, and the report says so explicitly. `[CODE]`

## 6.2 Cytoplasmic inheritance and linkage → one locus, not a genome scan

Plant chloroplasts and mitochondria are commonly inherited uniparentally, often
maternally, but paternal leakage and biparental inheritance occur. This repository
does not establish the inheritance mode in *Dudleya*. For these analyses, each
organelle should be treated as a separate, largely linked **cytoplasmic locus**,
not as a set of independent nuclear markers. Plant mitochondrial genomes can
recombine among repeats (§6.4), so “largely linked locus” is more accurate than
“strictly non-recombining molecule.” See reviews of [plant organelle
inheritance](https://www.nature.com/articles/hdy199419) and [plant mitochondrial
recombination](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/j.1469-8137.2005.01492.x).

The computational consequence is large and easy to forget: **the entire cpDNA
alignment is effectively one marker, and the entire mtDNA alignment is one
marker.** The pipeline computes hundreds of cpDNA SNPs and a tree with 275 tips,
but those SNPs are not independent — they mostly share one organelle genealogy.
So:

- A cpDNA tree is a *gene tree* of one cytoplasmic locus, not a species tree.
  High bootstrap support means the data strongly support that
  *organelle-haplotype* topology, not that it is the species history.
- Statistical "sample sizes" are misleading. 2,015 cpDNA SNPs are not 2,015
  independent observations of population history; they are one genealogy sampled
  at 2,015 points.
- cpDNA and mtDNA are two separate cytoplasmic markers, so comparing their trees
  is meaningful. They are not automatically statistically independent: shared
  transmission can couple their histories. Agreement is concordance between the
  markers; disagreement is discordance to investigate, neither of which by itself
  establishes nuclear history.

This single fact — organelle = one linked locus — is the reason the book's
recurring caveat exists. It is developed fully in [Chapter
17](./17-uncertainty-bias-and-limits.md). `[BIO]`

## 6.3 The chloroplast inverted repeat → the duplicate-copy mask

Most land-plant chloroplast genomes contain a large **inverted repeat (IR)**: two
near-identical copies of a multi-kilobase region in opposite orientation,
separated by a large and a small single-copy region. In this project's
normalized chloroplast reference the two IR copies sit at roughly
`82091–107826` and `124539–150274`. `[RESULT]`

Two identical copies are poison for two different computations:

1. **Mapping.** A read from the IR could belong to either copy, so it maps
   ambiguously (low MAPQ) or arbitrarily to one copy.
2. **Variant counting.** If you keep both near-identical copies, homologous IR
   signal can be duplicated and strongly correlated — inflating SNP counts and
   distorting distance or diversity statistics.

The pipeline's answer is to **keep exactly one IR copy** for population genetics.
Stage 05 builds `cpdna_population_sites.bed` as the whole chloroplast *minus* the
later IR copy, so the population-genetic alignment is 124,538 bp rather than the
full 150,274 bp ([Chapter 9](./09-masks-alignment-and-sample-qc.md)). Coverage QC
still uses the full reference, because for "did this sample sequence well?" you
want to see the whole molecule. `[CODE]` This is the concrete meaning of the
book's rule that *QC tracks and population-genetic tracks are not
interchangeable*.

## 6.4 Mitochondrial repeats and rearrangement → the high-confidence unique track

Plant mitochondrial genomes are the awkward organelle: large, repeat-rich, prone
to recombination among repeats, and often carrying sequence of plastid or
nuclear origin. Reads from repetitive or paralogous mtDNA regions map with low
confidence or to the wrong place, which would create false variants.

The pilot investigation (Stage 03) quantified this by comparing coverage at
permissive MAPQ against coverage at high MAPQ: the high-MAPQ, unique-placement
breadth was much lower than the permissive-MAPQ depth suggested. So the pipeline
does **not** trust the whole mitochondrial reference for variant calling. It
restricts mtDNA variants and population genetics to
`mtdna_high_confidence_unique_regions.bed` — currently just two high-MAPQ
consensus intervals totaling **44,930 bp**, a fraction of the 243,359 bp
reference. `[RESULT]`

The downstream price is visible in the numbers: mtDNA yields only **146** filtered
SNPs versus cpDNA's **2,015** ([Chapter 1, §1.5](./01-data-flow-map.md)). Restricting
analysis to 44,930 uniquely mappable bases contributes directly to that smaller
marker set, but the count alone cannot establish how variable the two organelles
are per callable base. Fewer sites give the mtDNA analyses less site-level
resolution, and the book flags that wherever mtDNA appears. `[BIO]`

## 6.5 Callability → depth, missingness, and `N`

A base is only usable if enough reads covered it confidently. This pipeline
treats *callability* as a first-class idea at three levels:

- **Per base, per sample.** In the callable-site consensus (Stage 11), any base
  below the minimum depth is written `N`. So the alignment records not just
  *what* base a sample has, but *whether we could see it at all*.
- **Per site, across samples.** Variant filtering drops any site genotyped in
  fewer than 80% of samples (`F_MISSING<=0.2`), because a variant seen in only a
  handful of samples is more likely an artifact than signal.
- **Per sample, overall.** Sample QC (Stage 06) excludes samples whose organelle
  breadth is too low; three low-input samples were dropped this way, leaving 275.

Across all 275 samples, the mtDNA callable consensus has 31,313 missing cells
(**0.2534%** of 275 × 44,930) and cpDNA has 127,485 (**0.3722%** of
275 × 124,538). `[RESULT]` Missingness is low overall but may be concentrated in
particular samples or sites. Every downstream method still has to handle it: PCA
mean-imputes missing sites, Fst skips sites where a population has no data, and
ADMIXTURE encodes missing as `0 0`. `[BIO]`

## 6.6 Reference bias → what mapping to one genome can and cannot see

Every sample is mapped to a *single* Dudleya organelle reference. Reads too
divergent from that reference map poorly or not at all, so the analysis is
blind to variation in regions that differ a lot from the reference individual —
this is **reference bias**. For organelle work it is usually mild (organelle
genomes are conserved within a genus) but it is not zero, and it interacts with
the mtDNA repeat problem: divergent *and* repetitive regions are doubly hard.

The pipeline's honesty about this is structural: it restricts calling to regions
where mapping is trustworthy (§6.3–6.4), records exactly which regions those are
as BED files, and reports breadth so you can see where coverage thins. It cannot
*remove* reference bias with one reference, and the book does not pretend
otherwise. `[BIO]`

## 6.7 Why cpDNA and mtDNA share a reference but never share an analysis

The reference FASTA has both organelles as two records so that one mapping pass
sorts each read to its correct molecule (via `samtools idxstats` per record).
But from that point on they are two experiments:

- separate BED tracks, separate VCFs, separate alignments;
- separate trees, PCA, clustering, and Fst;
- separate best-K and separate SNP counts.

They are combined only in the final *narrative* report, and even there they are
compared, not merged. The biological justification is everything in this chapter:
they differ in size, repeat content, mutation rate, and callable fraction, so
pooling them would blend two very different data qualities. Keeping them apart is
also what makes cpDNA-vs-mtDNA discordance detectable, which is one of the more
interesting things two organelle markers can tell you. `[BIO]`

## 6.8 The biology-to-code map

Keep this table nearby as you read Part 2. Each biological fact on the left is a
concrete decision on the right.

| Biology | Pipeline decision | Where |
|---|---|---|
| Organelles are haploid | `bcftools call --ploidy 1`; single-allele genotypes | Stage 08 |
| ADMIXTURE assumes diploid | pseudo-diploid homozygous encoding | Stage 16/18 |
| cpDNA inverted repeat duplicates sequence | keep one IR copy: `cpdna_population_sites.bed` (124,538 bp) | Stage 05 |
| mtDNA repeats map ambiguously | restrict to `mtdna_high_confidence_unique` (44,930 bp) | Stage 03 → 05 |
| Coverage QC ≠ variant trust | separate full-reference QC track vs population track | Stage 05 |
| A base needs enough reads | min-depth mask → `N` in consensus | Stage 11 |
| A site needs enough samples | `F_MISSING<=0.2` filter | Stage 09 |
| Each organelle = one largely linked cytoplasmic locus | trees are gene trees; no nuclear/species claim | Part 3 |

You now have the biology, the Python, the patterns, the tools, and the formats.
Part 2 walks the pipeline stage by stage, and every chapter there will reach back
to something in Part 1.

> Next: [Chapter 7 — Manifest and Reference Preflight (Stages 00–01)](./07-manifest-and-reference-preflight.md)
