# Chapter 12 — PCA and Admixture-Style Clustering (Stages 15, 16, 18)

> Part 2 of 4 · Pipeline Walkthrough · Prev:
> [Phylogenetic Trees](./11-phylogenetic-trees.md) · Next: [Population
> Fst](./13-population-fst.md)

Both stages consume the compact SNP-only alignment from Stage 10 and summarize
population structure — PCA as a scatter of samples, ADMIXTURE as stacked ancestry
bars. Building them is here; reading them is [Chapter
16](./16-pca-clustering-fst-interpretation.md).

## 12.1 The question

*Without assuming any grouping, do samples cluster by species and population, and
how many organelle-haplotype clusters best describe the data?*

## 12.2 Stage 15 — PCA

[`pca_analysis.py`](../pca_analysis.py) (runner:
[`../scripts/run_pca_analysis.py`](../scripts/run_pca_analysis.py)) turns the
SNP FASTA into a numeric matrix, runs a 2-component PCA, and writes coordinates,
variance, and figures for each organelle.

### Encoding SNPs as numbers

`build_haploid_snp_matrix` walks each SNP column, keeps only columns with at
least two observed alleles, integer-codes the alleles, and **mean-imputes**
missing values:

```python
alleles = sorted({base for base in site_bases if base in BASES})
if len(alleles) < 2:
    continue                                 # monomorphic column: drop
allele_codes = {base: float(index) for index, base in enumerate(alleles)}
encoded = np.array([allele_codes.get(base, np.nan) for base in site_bases])
site_mean = float(np.nanmean(encoded))
encoded[np.isnan(encoded)] = site_mean       # N -> column mean
columns.append(encoded)
```

So at a site with alleles `A` and `T`, `A→0.0`, `T→1.0`, and any `N` becomes the
column's mean (a value between 0 and 1). A test on a three-sample FASTA
`AAN / ATG / TTT` asserts the imputed cell equals `0.5`. `[TEST]` Mean imputation
is a pragmatic choice — it lets PCA run without dropping every partially-missing
site — but it does pull imputed samples toward the center, which matters for
interpretation ([Chapter 16](./16-pca-clustering-fst-interpretation.md)).

### The PCA itself

`run_pca` mean-centers the matrix and calls scikit-learn:

```python
centered = matrix - matrix.mean(axis=0)
pca = PCA(n_components=2)
coordinates = pca.fit_transform(centered)
return coordinates, pca.explained_variance_ratio_
```

The observed variance explained: cpDNA **PC1 36.62%, PC2 14.65%**; mtDNA
**PC1 34.48%, PC2 14.06%**. `[RESULT]` Points are colored by `choose_plot_group`,
which prefers `species_popcode`, then popcode, then species, then naming profile —
so samples without resolved metadata still get a sensible label. The stage writes
coordinates and variance TSVs and PNG/PDF/SVG scatterplots. A test runs the whole
`run_one_pca` on a four-sample fixture and asserts all outputs exist. `[TEST]`

## 12.3 Stage 16/18 — ADMIXTURE-style clustering

[`admixture_analysis.py`](../admixture_analysis.py) (runner:
[`../scripts/run_admixture_analysis.py`](../scripts/run_admixture_analysis.py)).
Stage 16 is the single-run first pass; **Stage 18 is the final deliverable** with
five seeded replicates per K.

### Pseudo-diploid encoding

ADMIXTURE assumes diploids, so the haploid calls are duplicated into homozygous
pairs in a PLINK PED/MAP ([Chapter 5, §5.7](./05-bioinformatics-file-formats.md)).
Samples whose every SNP is missing are excluded and recorded:

```python
if all(base not in BASES for base in sequence):
    excluded_sample_ids.append(sample_id)
    continue
...
for base in sequence:
    allele = base if base in BASES else "0"
    genotype_fields.extend([allele, allele])       # A -> "A A", N -> "0 0"
```

Tests assert `TG → "T T G G"`, `AN → "A A 0 0"`, and that an all-missing sample
is dropped into the excluded table with reason `all_snp_genotypes_missing`.
`[TEST]`

### PLINK, the K sweep, and CV

`plink --make-bed` converts PED/MAP to the binary triple, then ADMIXTURE runs for
each K with cross-validation. The K ceiling is clamped so the run is always
valid:

```python
k_max_for_input = min(max_k, len(included_sample_ids) - 1, admixture_input.alignment_sites)
```

For each K and replicate, `run_admixture_for_k` runs `admixture --cv
--seed=<seed> -j4 <bed> <K>`, advancing the seed per replicate, and parses the CV
error from the log with a regex:

```python
match = re.search(r"CV error \(K=\d+\):\s*([0-9.eE+-]+)", log_text)
```

A test asserts the command is exactly `["admixture", "--cv", "--seed=42", "-j4",
"cpDNA.ped", "3"]` and that CV parsing reads `0.12345` from a log line. `[TEST]`

### Selecting K by mean CV error

`summarize_replicate_stability` groups replicates by (organelle, K), averages
their CV errors, and marks the K with the lowest *mean* as best:

```python
best_by_organelle[organelle] = min(organelle_items, key=lambda item: item[1])[0]
```

A test with two replicates each at K=2 (mean 0.31) and K=3 (mean 0.21) asserts
K=3 is `is_best_mean_k = yes` with mean `0.21000000`. `[TEST]` The best replicate
at the best K supplies the Q matrix for the structure plot.

### The result and its limits

Both organelles have their lowest tested mean CV error at **K = 8** in the
five-replicate Stage 18: cpDNA 0.08899 (SD 0.01449), mtDNA 0.12644 (SD 0.02207).
`[RESULT]` For both organelles the mean CV error decreases monotonically from K=1
to the upper boundary K=8. The tested range therefore did **not bracket an
interior optimum**; it does not show that eight is correct, or that at least eight
biological groups exist.

There is a deeper limitation: ADMIXTURE is a diploid population-model program,
while this stage duplicates haploid organelle calls into homozygous pseudo-diploid
genotypes and supplies many linked sites from each cytoplasmic locus. Its model
does not explicitly account for that linkage. Treat the K sweep and Q bars as an
exploratory similarity visualization under violated assumptions, not as validated
ancestry proportions, haplotype assignments, or a count of biological groups.
The [original ADMIXTURE paper](https://faculty.eeb.ucla.edu/Novembre/AlexanderEtAl_GR_2009.pdf)
explicitly notes that the model does not account for linkage disequilibrium. This
is developed in [Chapter 16](./16-pca-clustering-fst-interpretation.md).
`[BIO]`

The structure plot sorts samples by metadata group so related samples sit
together, and colors are the `tab20` colormap. Stage 18 also writes a CV-vs-K
plot per organelle.

## 12.4 The Python concepts here

- **NumPy arrays and `np.nanmean`** for column encoding and imputation.
- **scikit-learn `PCA`** — mean-center, `fit_transform`, `explained_variance_ratio_`.
- **Regex extraction** of a number from tool stdout.
- **`min(..., key=lambda ...)`** to select an extremum by a computed key.
- **Running a tool with `cwd=output_dir`** and renaming its generic outputs — a
  workaround for a tool that names files after its input.
- **`os.environ.setdefault("MPLCONFIGDIR", ...)`** for headless matplotlib.

## 12.5 Failure modes

- **Missing `admixture` or `plink`** → `AdmixtureAnalysisError`. `[CODE]`
- **No polymorphic SNP columns** → `PcaAnalysisError` (nothing to decompose).
  `[CODE]`
- **Fewer than 2 samples/sites** → PCA raises. `[CODE]`
- **All-missing sample** → excluded from ADMIXTURE, recorded, not fatal (unless
  *no* informative samples remain, which raises). `[CODE]`
- **ADMIXTURE writes no Q/P files** → `AdmixtureAnalysisError`. `[CODE]`
- **Statistical failure: over-reading structure.** PCA axes and cluster counts
  from 146 mtDNA SNPs are noisy; pseudo-diploid ADMIXTURE plots are *haplotype*
  clustering, not nuclear admixture. Both cautions are Chapter 16. `[BIO]`

## 12.6 Exercises

1. **Trace.** A SNP column across four samples is `A, A, T, N`. What allele codes
   does `build_haploid_snp_matrix` assign, and what value replaces the `N`?
2. **Predict.** A SNP column is `G, G, G, G` (monomorphic). Is it kept in the PCA
   matrix? What about a column `A, A, N, N`?
3. **Predict.** `build_admixture_command("admixture", Path("mtDNA.ped"), k=5,
   seed=100, threads=8)` returns what list?
4. **Modify.** You want to color the PCA by `naming_profile` only. Which one
   function do you change, and what would `choose_plot_group` return for a sample
   with species and popcode both blank?
5. **Debug.** ADMIXTURE has its lowest tested CV error at K=1. What CV-table shape
   would produce that, and what limited statement can you make under this model?
6. **Interpret.** Both organelles have their lowest tested CV error at K=8, the
   maximum tested. Explain why the monotonic curve fails to identify a biological
   group count, then name the ploidy and linkage assumptions that further limit
   interpretation.

Solutions in [Chapter 19](./19-solutions.md).

> Next: [Chapter 13 — Population Fst and Diversity (Stage 17)](./13-population-fst.md)
