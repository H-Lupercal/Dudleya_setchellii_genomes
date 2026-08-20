# Chapter 13 — Population Fst and Diversity (Stage 17)

> Part 2 of 4 · Pipeline Walkthrough · Prev: [PCA and
> Clustering](./12-pca-and-clustering.md) · Next: [The Tool
> Audit](./14-tool-audit.md)

This is the pipeline's one hand-written population-genetics stage: no external
popgen tool, just Python and `collections.Counter`. It computes pairwise Fst
between populations and per-population diversity from the SNP alignment. Reading
the numbers is [Chapter 16](./16-pca-clustering-fst-interpretation.md); here we
read the math.

## 13.1 The question

*Among the populations with known codes, how differentiated is each pair
(pairwise Fst), and how much variation does each population hold (haplotype and
nucleotide diversity, private variants)?*

## 13.2 The files

[`population_genetics.py`](../population_genetics.py) (runner:
[`../scripts/run_population_genetics.py`](../scripts/run_population_genetics.py))
reads the Stage 10 SNP FASTA and the sample metadata, and writes, per organelle:
`*.pairwise_fst.tsv`, `*.population_summary.tsv`, a stage summary, and a
`report.md`.

## 13.3 Only resolved populations

Fst needs named populations, so the stage keeps only samples whose `popcode` is
non-empty and groups them:

```python
def group_sequences_by_population(records, metadata):
    groups = defaultdict(list)
    for sample_id, sequence in records:
        popcode = population_code_for_sample(metadata.get(sample_id, {}))
        if popcode:
            groups[popcode].append((sample_id, sequence))
    return dict(sorted(groups.items()))
```

Samples from the initial DU-only batches have no popcode, so they contribute to
the alignments, PCA, and trees but are **omitted** from population summaries.
The observed run has **34** metadata-resolved populations, giving
34 × 33 / 2 = **561** pairwise comparisons per organelle. `[RESULT]` If fewer than
two populations resolve, the stage raises `PopulationGeneticsError` — you cannot
compute between-population differentiation with one group.

## 13.4 The Fst math

The pipeline uses a Nei-style gene-diversity Fst, averaged over informative
sites. Gene diversity at a site is `1 − Σ p²` (the chance two random alleles
differ):

```python
def gene_diversity(counts: Counter[str]) -> float:
    total = sum(counts.values())
    if total <= 1:
        return 0.0
    return 1.0 - sum((count / total) ** 2 for count in counts.values())
```

For a pair of populations, at each site the code computes each population's
within-diversity (`h1`, `h2`), the combined total diversity (`ht`), the average
within (`hs = (h1 + h2)/2`), and accumulates a numerator and denominator:

```python
numerator   += max(ht - hs, 0.0)
denominator += ht
...
return numerator / denominator, informative_sites
```

So Fst = Σ(H_T − H_S) / Σ H_T over informative sites — high when populations are
differentiated (their combined diversity exceeds their within-group diversity),
zero when they are identical. Two tests pin the behavior: two populations fixed
for different alleles (`AA,AA` vs `TT,TT`) give Fst **1.0** across 2 informative
sites, and haplotype diversity of `AA, AA, AT` is exactly **2/3**. `[TEST]`
Sites where either population has no data, or where the combined site is
monomorphic, are skipped (`if not counts1 or not counts2: continue`) — missingness
is handled by exclusion, not imputation, unlike PCA.

## 13.5 The diversity statistics

Each population also gets four summary numbers:

- **Haplotype count** — `len(set(sequences))`, the number of distinct whole
  sequences (valid because each sample is one haplotype, [Chapter
  6](./06-organelle-biology.md)). Important implementation caveat: this exact
  string comparison includes `N`, so two otherwise identical records with
  different missing-data patterns count as different haplotypes.
- **Haplotype diversity** — `(n/(n−1)) · (1 − Σ f²)`, the probability two random
  samples are different sequence strings, bias-corrected for sample size. Because
  `N` is retained in those strings, missing-data patterns can inflate the count
  and diversity; inspect them before giving the statistic a biological reading.
- **Nucleotide diversity** — the average per-site difference over all sample
  pairs, counting only positions where both samples have a real base:

  ```python
  for base1, base2 in zip(seq1, seq2):
      if base1 not in BASES or base2 not in BASES:
          continue
      compared_sites += 1
      if base1 != base2:
          differences += 1
  ```

- **Private variant sites** — sites where the population carries an allele seen
  in no other population (`pop_alleles - other_alleles`, set difference again).

These are computed with `itertools.combinations` over sample pairs and
`collections.Counter` over alleles — standard library only, no external popgen
dependency. That keeps the math auditable line by line, which is the whole point
of doing it in Python here.

## 13.6 The Python concepts here

- **`collections.Counter`** for allele frequencies and its `+` merge.
- **`itertools.combinations`** for all population pairs and all sample pairs.
- **Guarded division** everywhere (`if total <= 1`, `if denominator <= 0`) so
  degenerate inputs return 0.0 instead of crashing.
- **Set difference** for private variants.
- **`max(ht - hs, 0.0)`** to floor a differentiation term at zero.

## 13.7 The result, stated honestly

The stage produces 561 pairwise Fst values and 34 population diversity rows per
organelle. `[RESULT]` What Fst here **is**: a relative measure of organelle
haplotype differentiation between two populations, on the trusted sites, using a
transparent Nei-style estimator.

What it is **not**: a nuclear Fst, and not comparable to Fst values from other
studies computed with different estimators or on nuclear SNPs. It is also
sensitive to sample size and to the very different SNP counts of the two
organelles (2,015 cpDNA vs 146 mtDNA) — an mtDNA Fst rests on far fewer
informative sites. How to read these numbers, and what over-reading them looks
like, is [Chapter 16](./16-pca-clustering-fst-interpretation.md). `[BIO]`

## 13.8 Failure modes

- **Fewer than two resolved populations** → `PopulationGeneticsError`. `[CODE]`
- **Missing SNP FASTA** → raised by `read_population_inputs`. `[CODE]`
- **A population with a single sample** → diversity statistics return 0.0 by
  construction (`if n <= 1`), not an error — but a single-sample "population" is
  a weak basis for any claim. `[BIO]`
- **Limited mtDNA site count** → the calculation is mathematically valid, but 146
  available SNPs provide less site-level resolution than 2,015. Judge each pair
  using its actual `informative_sites` count rather than dismissing or accepting
  mtDNA Fst from the headline count alone. `[BIO]`

## 13.9 Exercises

1. **Trace.** Population A has sequences `AA, AA`; population B has `TT, TT`. Work
   `gene_diversity` for A, B, and the combined set at site 1, then the resulting
   Fst. Confirm it matches the test's `1.0`.
2. **Predict.** Population C is `AA, AT, TT`. What is its haplotype count and
   `compute_haplotype_diversity`?
3. **Predict.** With 34 populations, how many rows are in the pairwise Fst table?
   With 10 populations?
4. **Modify.** You want to exclude populations with fewer than 3 samples from Fst.
   Where in `run_one_population_summary` would you filter `groups`, and what
   effect would it have on the comparison count?
5. **Debug.** An mtDNA pairwise Fst is exactly 0.0 for two populations you expect
   to differ. Name two reasons (missingness, monomorphic sites) and the column in
   the pairwise table that would help you tell which.
6. **Interpret.** cpDNA Fst between two populations is 0.4; mtDNA Fst between the
   same pair is 0.0. Give a data-quality explanation and a biological explanation,
   and say which you would rule out first and how.

Solutions in [Chapter 19](./19-solutions.md).

> Next: [Chapter 14 — The Bioinformatics Tool Audit (Stage 13)](./14-tool-audit.md)
