# Chapter 10 — From Reads to Alignments (Stages 08–11)

> Part 2 of 4 · Pipeline Walkthrough · Prev: [Masks,
> Alignment, and Sample QC](./09-masks-alignment-and-sample-qc.md) · Next:
> [Phylogenetic Trees](./11-phylogenetic-trees.md)

This is the pipeline's core. Four stages turn 275 BAMs into the two alignment
products everything downstream consumes: a compact SNP-only alignment (for PCA,
clustering, Fst) and a full callable-site consensus (for trees). Read this
chapter carefully; the interpretation chapters lean on it.

## 10.1 The question

*Where do the samples differ from each other and from the reference, which of
those differences do we trust, and how do we write them as alignments that also
record where we could not see?*

## 10.2 The four stages at a glance

| Stage | Module | In | Out (cpDNA / mtDNA) |
|---|---|---|---|
| 08 variant calling | [`variant_calling.py`](../variant_calling.py) | BAMs + tracks | raw VCF; 2,475 / 190 records |
| 09 variant filtering | [`variant_filtering.py`](../variant_filtering.py) | raw VCF | filtered SNP VCF; 2,015 / 146 |
| 10 SNP alignment | [`snp_alignment.py`](../snp_alignment.py) | filtered VCF | SNP-only FASTA; 2,015 / 146 sites |
| 11 callable consensus | [`callable_consensus.py`](../callable_consensus.py) | raw+filtered VCF, BED, depth | full FASTA; 124,538 / 44,930 sites |

`[RESULT]` Notice the two alignments have very different shapes: the SNP
alignment is *only the variable sites*, while the callable consensus is *every
trusted position*, most of them invariant or `N`. They exist for different
downstream methods.

## 10.3 Stage 08 — haploid variant calling

Run by [`../scripts/run_variant_calling.py`](../scripts/run_variant_calling.py).
The stage reads `included_samples.tsv`, resolves each sample's BAM by its
filesystem-safe name, and restricts calling to the population track for each
organelle:

```python
wanted = {"cpdna_population_sites": "cpDNA",
          "mtdna_high_confidence_unique": "mtDNA"}
```

Only those two tracks are used — a test asserts the QC tracks are ignored here.
`[TEST]` For each organelle it runs the `bcftools mpileup | call` pipe with
**`--ploidy 1`**, then `reheader` (to restore real sample names) and `index`
([Chapter 4, §4.4](./04-shell-and-external-tools.md)). Calling defaults are
`min_mapq = 20`, `min_baseq = 20`, `max_depth = 10000` — deliberately stricter
than the mapping stage so variants come only from confidently placed, confidently
read bases. `[CODE]`

The observed raw call counts are **2,475** cpDNA and **190** mtDNA variant
records across 275 samples. `[RESULT]` The gap already tells the organelle-size
story from [Chapter 6](./06-organelle-biology.md): a large, fully-used chloroplast
track versus a small, trusted mtDNA track.

### Safe filenames and sample identity

Sample IDs like `KEEP/ONE` contain characters illegal in filenames, so BAMs are
stored under `safe_sample_name` (slashes → underscores) while the *real* ID is
preserved for the VCF header via `reheader -N`. A test builds a `KEEP/ONE` sample
and asserts `safe_sample_id == "KEEP_ONE"` and the BAM path uses the safe name.
`[TEST]` This split — safe name on disk, real name in the data — recurs in the
consensus stage too.

### The smoke run

Stage 08 is the one stage with a `smoke` run label: five named samples, run first
to validate the calling machinery before spending compute on all 275
([Chapter 1, §1.4](./01-data-flow-map.md)). The smoke outputs live beside the
primary ones but never feed the deliverables.

## 10.4 Stage 09 — filtering to trustworthy SNPs

Run by
[`../scripts/run_variant_filtering.py`](../scripts/run_variant_filtering.py). It
reads the Stage 08 summary to find the raw VCFs and applies one `bcftools view`
filter ([Chapter 4, §4.5](./04-shell-and-external-tools.md)):

```text
-m2 -M2 -v snps  --min-ac 2:minor  -i 'F_MISSING<=0.2'
```

Biallelic, SNP-only, minor allele seen at least twice, genotyped in ≥80% of
samples. The four thresholds are `DEFAULT_MAX_MISSING_FRACTION = 0.2` and
`DEFAULT_MIN_MINOR_ALLELE_COUNT = 2`. The effect: **2,475 → 2,015** cpDNA and
**190 → 146** mtDNA. `[RESULT]` About a fifth of raw cpDNA calls and a quarter of
raw mtDNA calls are dropped as multiallelic, indel, singleton, or too-missing —
the difference between "something varied here" and "a variant we would defend."

`FilterInput` derives the filtered path from the raw path by string replacement
(`.raw.vcf.gz` → `.filtered.vcf.gz`), and a test pins that derivation. `[TEST]`

## 10.5 Stage 10 — the SNP-only alignment

Run by [`../scripts/run_snp_alignment.py`](../scripts/run_snp_alignment.py). This
stage is pure Python: it reads the filtered VCF and, for each variant site,
writes one base per sample using `genotype_to_base` ([Chapter 5,
§5.6](./05-bioinformatics-file-formats.md)):

```python
for sample, genotype_field in zip(sample_names, columns[9:], strict=True):
    sequence_parts[sample].append(genotype_to_base(genotype_field, ref, alt))
```

The output is a FASTA where each record is one sample and each column is one SNP:
2,015 columns for cpDNA, 146 for mtDNA. Missing genotypes become `N`. A test
feeds a two-site, three-sample VCF and asserts the exact sequences `"AC"`,
`"GT"`, `"NN"`. `[TEST]` The stage also enforces that every site is truly
biallelic (single-character REF and ALT), raising `SnpAlignmentError` otherwise —
a belt-and-suspenders check on top of the Stage 09 filter. This compact alignment
is the input for PCA, clustering, and Fst.

## 10.6 Stage 11 — the full callable-site consensus

Run by
[`../scripts/run_callable_consensus.py`](../scripts/run_callable_consensus.py).
This is the most intricate stage in the pipeline, and worth slowing down for. Its
goal: a FASTA covering *every position in the population track*, where each
sample's base is the reference by default, the called allele where one passed
filtering, and `N` wherever the site failed filtering or the sample lacked depth.

The algorithm, from `build_callable_consensus`, in five moves:

1. **Template.** Build the per-position reference sequence over the track's BED
   intervals, and a `coordinate_to_index` map from `(chrom, position)` to a
   column index. Every sample starts as a copy of this reference template.

   ```python
   for position in interval.positions_1based:
       base = record[position - 1].upper()
       coordinate_to_index[(interval.chrom, position)] = len(bases)
       bases.append(base)
   ```

2. **Failed-site set.** Compute which raw-variant sites did *not* survive
   filtering — set difference again ([Chapter 2, §2.8](./02-python-essentials.md)):

   ```python
   failed_site_indexes = {coordinate_to_index[k]
                          for k in raw_site_keys - filtered_site_keys}
   ```

   These are positions where the caller saw *something* but the filter rejected
   it; they will be masked to `N` because we neither trust the variant nor the
   reference there.

3. **Overlay filtered SNPs.** At each filtered biallelic SNP, replace every
   sample's base with its `genotype_to_base` call.

4. **Depth mask.** For each sample, read its Stage 06 depth file, and write `N` at
   every failed-site index and at every position not covered to at least
   `min_depth` (default 1):

   ```python
   for index in failed_site_indexes:      sequence[index] = "N"
   for index in range(len(sequence)):
       if index not in covered:           sequence[index] = "N"
   ```

5. **Assemble** the per-sample strings and count missing bases.

The observed products: cpDNA is 275 × **124,538** with **459** masked
failed-variant sites and **127,485** total missing cells; mtDNA is 275 ×
**44,930** with **44** masked failed sites and **31,313** missing cells.
`[RESULT]` The denominator is the full sample-by-site matrix, not alignment width:
cpDNA missingness is **0.3722%** and mtDNA missingness is **0.2534%**. These low
overall rates can still conceal samples or sites with locally poor callability,
which is why downstream methods must treat `N` explicitly ([Chapter 6,
§6.5](./06-organelle-biology.md)).

### A subtlety worth naming

The depth files that define callability come from Stage 06, produced at the
mapping default `min_mapq = 0` (permissive). But the *variants* overlaid here
come from Stage 08 at `min_mapq = 20` (strict). So a consensus position can be
"callable" (had ≥1× permissive-MAPQ coverage) while carrying only the reference
base because no high-MAPQ variant was called there. This is a reasonable design —
be generous about "we saw this position," strict about "we call a variant here" —
but it is a real asymmetry, and you should know it when interpreting the trees
built from this alignment. `[CODE]`

The whole Stage 11 algorithm is traced by hand on a six-base example in the
[capstone, Chapter 23](./23-capstone-sample-trace.md); the test that pins it
(`"CTTNNG"`/`"CGTANG"`) is the same example.

## 10.7 The Python concepts here

- **`bcftools` command construction** as argument lists, pinned by
  command-builder tests ([Chapter 3, §3.8](./03-reusable-code-patterns.md)).
- **Set difference** to compute failed-variant sites.
- **Coordinate→index maps** to align VCF positions, BED positions, and alignment
  columns.
- **The shared `genotype_to_base`** used identically in Stages 10 and 11.
- **Filesystem-safe naming** decoupled from biological sample identity.
- **String-replacement path derivation** (`.raw` → `.filtered`), guarded by a
  test.

## 10.8 Failure modes

- **Missing BAM/index for an included sample** → `VariantCallingError`. `[CODE]`
- **Requested sample absent from the included set** → `VariantCallingError`
  listing the missing IDs. `[CODE]`
- **Missing population track** → `VariantCallingError`. `[CODE]`
- **Raw or filtered VCF sample order ≠ included samples** →
  `CallableConsensusError` — a guard that the VCFs and the sample table agree.
  `[CODE]`
- **Missing depth file for a sample** → `CallableConsensusError`. `[CODE]`
- **A non-biallelic site in a "filtered" VCF** → `SnpAlignmentError`. `[CODE]`
- **Biological caution: `N` is not `A`.** Every downstream method must treat `N`
  as absence of data, not a fifth base. Overall missingness is low in both
  organelles, but unusually missing samples or sites can still distort methods
  that impute or compare incomplete observations. `[BIO]`

## 10.9 Exercises

1. **Trace.** A filtered VCF row is `chloroplast 20 . T C . PASS . GT 1 0 .`.
   What three bases does `genotype_to_base` produce for the three samples?
2. **Predict.** Raw calls include sites at positions 3 and 6; the filter keeps
   only position 3. In the callable consensus, what happens to every sample's
   base at position 6, and why?
3. **Predict.** You raise `--min-depth` from 1 to 5 in Stage 11. Directionally,
   how do `missing_bases` and the number of informative alignment columns change?
4. **Modify.** You want Stage 09 to allow 30% missingness. Which constant and
   which `bcftools view` flag change, and would you expect cpDNA filtered SNP
   count to go up or down?
5. **Debug.** Stage 11 raises "Raw VCF sample order does not match included
   samples." Give two plausible causes and the files you would compare.
6. **Interpret.** cpDNA has 2,015 SNPs in a 124,538-bp alignment; mtDNA has 146
   in 44,930 bp. Compute the rough SNP density for each. What does the difference
   imply for the resolution of the mtDNA tree versus the cpDNA tree?

Solutions in [Chapter 19](./19-solutions.md).

> Next: [Chapter 11 — Phylogenetic Trees (Stages 12, 14, 19, 20)](./11-phylogenetic-trees.md)
