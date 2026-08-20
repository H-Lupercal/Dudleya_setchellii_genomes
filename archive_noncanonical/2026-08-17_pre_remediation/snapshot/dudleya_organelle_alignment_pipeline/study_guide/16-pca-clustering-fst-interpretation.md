# Chapter 16 — Reading PCA, Clustering, and Fst

> Part 3 of 4 · Interpretation · Prev: [Reading the
> Trees](./15-phylogenetics-interpretation.md) · Next: [Uncertainty, Bias, and
> Limits](./17-uncertainty-bias-and-limits.md)

Three summaries of the same SNP alignments, three ways to over-read them. This
chapter is about extracting what each one honestly supports.

## 16.1 PCA: relative, unsupervised, mean-imputed

The PCA scatter places each sample by its first two principal components of the
SNP matrix. Read it as *"which samples resemble which,"* with four constraints:

- **Axes are relative and unitless.** PC1 is just the direction of greatest
  variance. "Far apart on PC1" means "differ along the dominant axis of organelle
  SNP variation," nothing more. The variance captured is modest: cpDNA
  PC1 36.62% / PC2 14.65%, mtDNA PC1 34.48% / PC2 14.06% — so roughly half the
  variation lives in components you are not looking at. `[RESULT]`
- **It is unsupervised.** The colors come from metadata
  (`choose_plot_group`, [Chapter 12](./12-pca-and-clustering.md)), but the
  *positions* do not. When colored groups separate cleanly, that is a real
  finding: the organelle SNPs alone recover the labeling. When they overlap, the
  organelle data simply do not distinguish those groups. `[BIO]`
- **Missing data is mean-imputed.** Samples with many `N`s get pulled toward the
  center of each imputed axis, so a sample near the origin may be *uninformative*
  (heavily missing) rather than *intermediate* (genuinely admixed-looking). Check
  a central sample's callable fraction before interpreting its position. `[BIO]`
- **Sign and rotation are arbitrary.** PCA can flip an axis run to run; "left vs
  right" carries no meaning, only "same side vs opposite side."

The mtDNA PCA, built from 146 SNPs, is far coarser than the cpDNA PCA from 2,015.
Fewer sites means fewer resolvable directions, so read the mtDNA scatter as a
low-resolution sketch.

## 16.2 ADMIXTURE: an exploratory, assumption-limited projection

This is the most mislabel-prone plot in the pipeline, and the report says so in
its own text. Two things to internalize before you look at a single bar:

**It is not nuclear admixture, and its organelle interpretation is exploratory.**
The haploid calls were duplicated into pseudo-diploid homozygotes to satisfy
ADMIXTURE's diploid input model ([Chapter 5,
§5.7](./05-bioinformatics-file-formats.md)). Sites within each organelle are also
strongly linked, whereas the model does not explicitly account for linkage. A
“50% cluster A, 50% cluster B” bar is therefore neither evidence of a hybrid nor
a validated biological haplotype-assignment proportion. At most, it visualizes
how the misspecified model represents that sample relative to fitted components.
`[BIO]`

**K = 8 is a boundary result, not a lower bound on group count.** Both organelles
have their lowest tested mean CV error at K=8, the largest tested K, and CV error
decreases monotonically across K=1–8 ([Chapter 12,
§12.3](./12-pca-and-clustering.md)). `[RESULT]` The sweep therefore did not
bracket an interior optimum. It cannot distinguish among an optimum beyond the
range, continued overfitting, linkage effects, or other model mismatch, and it
does not establish “at least eight” biological groups. The replicate SDs (cpDNA
0.0145, mtDNA 0.0221 at K=8) describe variation among runs, not biological proof
of K.

**Reading the bars.** Samples are sorted by metadata group, so visual alignment
between colors and metadata can be described as an exploratory pattern. Do not
translate component colors into populations or ancestry. Ragged bars can reflect
limited signal, linked markers, model mismatch, or instability—not biological
mixture by default.

## 16.3 Fst: relative differentiation on trusted sites

Pairwise Fst ([Chapter 13](./13-population-fst.md)) is a relative measure of
organelle-haplotype differentiation between two populations. Read it comparatively
within this dataset, not as an absolute or cross-study number:

- **Higher means more differentiated** for that organelle, on the trusted sites,
  under this Nei-style estimator. A value near 1 means the populations are near
  fixed for different haplotypes; near 0 means they share the same variation.
- **Not comparable across studies or estimators.** Different Fst estimators and
  different marker sets give different numbers; do not line these up against a
  nuclear-SNP Fst from elsewhere. `[BIO]`
- **mtDNA Fst has fewer sites.** With only 146 SNPs, an mtDNA pairwise Fst may
  rest on fewer informative sites and be less precise. Overall missingness is
  low, so inspect each pair's `informative_sites` rather than assuming it is
  missingness-limited. cpDNA has more site-level resolution, but is a different
  linked locus and is not automatically biologically truer. `[BIO]`
- **Only 34 populations.** Samples without resolved population codes are excluded
  from Fst entirely, so the table describes the metadata-resolved subset, not all
  275 samples. `[RESULT]`

A useful cross-check: cpDNA Fst, cpDNA PCA separation, cpDNA tree clades, and
cpDNA exploratory ADMIXTURE components may *tell a consistent story* about the same
populations, because they are four views of the same 2,015-SNP alignment. Where
they disagree, suspect a method artifact (imputation, K ceiling, small
denominators) before a biological surprise.

## 16.4 The one-marker reminder, applied

PCA, clustering, and Fst here all describe *organelle* variation, with each
organelle treated as one largely linked cytoplasmic locus. Concordant cpDNA and
mtDNA patterns are worth reporting, but shared inheritance means they are not
automatically independent corroboration and they say nothing directly about
nuclear ancestry.
None of these plots can establish hybridization, gene flow in the nuclear genome,
or a species boundary on their own. That full argument is [Chapter
17](./17-uncertainty-bias-and-limits.md). `[BIO]`

## 16.5 What you may and may not say

| You may say | You may **not** say |
|---|---|
| "cpDNA PCA separates species X from Y along PC1." | "PC1 measures how much gene flow occurred." |
| "ADMIXTURE's lowest tested cpDNA CV error occurs at the boundary K=8." | "The data contain at least/exactly 8 populations or haplotype groups." |
| "Sample Z's exploratory Q bar contains multiple fitted components." | "Sample Z is an admixed hybrid or has validated mixed assignment." |
| "cpDNA Fst between A and B is high relative to other pairs here." | "The Fst of 0.4 is high by an absolute standard / matches nuclear Fst." |
| "mtDNA structure is low-resolution (146 SNPs)." | "mtDNA shows no structure, therefore none exists." |

## 16.6 Reading checklist

- [ ] For PCA: is a central sample intermediate, or just heavily missing?
- [ ] For ADMIXTURE: did you label the pseudo-diploid and linkage assumptions?
- [ ] For ADMIXTURE: did you report K=8 as a boundary result, not a group count?
- [ ] For Fst: are you comparing within this dataset and checking the number of
      informative sites?
- [ ] For all three: is your claim about *organelle* variation, tagged as such?

> Next: [Chapter 17 — Uncertainty, Bias, and the Limits of These Claims](./17-uncertainty-bias-and-limits.md)
