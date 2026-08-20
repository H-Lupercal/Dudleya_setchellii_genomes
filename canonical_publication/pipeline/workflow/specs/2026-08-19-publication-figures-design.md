# Canonical Publication Figures Design

Status: approved by the user on 2026-08-19.

## Purpose

Render a reproducible, publication-ready figure suite exclusively from fingerprint-validated canonical outputs. Figures are explanatory views of tabular/tree results, never independent analysis products, and no renderer may read from `archive_noncanonical/`.

## Architecture and data flow

`run_pca.py` and `run_haplotypes.py` produce only numerical/tabular analysis outputs. A dedicated `render_figures.py` stage runs after ADMIXTURE and consumes validated states from references, QC, PCA, haplotypes, population genetics, trees, and supplementary ADMIXTURE. It writes PNG, PDF, and SVG versions plus a figure manifest beneath `canonical_publication/reports/figures/<run-id>/` and records one checksummed provenance state.

The runner order is:

```text
... -> trees -> treecheck -> admixture -> figures -> reports
```

Reports fingerprint the figure state. Resume is accepted only when the figure fingerprint, upstream fingerprints, and output checksums agree.

## Visual system

- Use a fixed Okabe-Ito-derived, colorblind-safe five-taxon palette throughout PCA, haplotypes, and trees.
- Encode population with explicit text labels or annotation strips, not dozens of difficult-to-distinguish colors.
- Use a signed diverging FST scale centered at zero, with missing estimates in neutral gray.
- Use separate cluster colors for supplementary ADMIXTURE so they cannot be mistaken for taxon identity.
- Use deterministic layouts and record the layout seed or deterministic algorithm in provenance.
- Titles, captions, and legends state organelle, filtering basis, and supplementary status where applicable.
- Figures must remain legible in grayscale through labels, borders, and panel structure.

## Figure suite

1. **Reference and callability maps.** Chloroplast is circular because the selected molecule is circularly represented; mitochondria is linear because circularity is not supported. Annotation and mask tracks are labeled, and the mitochondrial panel distinguishes repeat exclusion from read-backed high-confidence sequence.
2. **Preprocessing and QC.** Show read retention and organelle-specific DP5 breadth/callability. Eligibility is visibly tied to the configured 80% DP5 threshold.
3. **Organelle-specific PCA.** Color samples by taxon, label population centroids with population codes, show explained variance, and state that only MAC>=2 markers were used.
4. **Organelle-specific haplotype networks.** Use exact mutational distances to calculate a deterministic network layout, label edges with mutational distance, scale nodes by sample count, and draw taxon-composition pies. Population detail remains in canonical assignment tables.
5. **Signed Hudson FST heatmaps.** Display all eligible population pairs, preserve negative estimates, center the scale at zero, and show unavailable values in gray.
6. **Primary unrooted organelle trees.** Render chloroplast and mitochondrial trees separately with branch lengths, support labels, and taxon-colored tips. Do not imply a root. The concatenated tree remains supplementary and is not promoted into a primary panel.
7. **Supplementary ADMIXTURE sensitivity figures.** Show replicate CV distributions/means and the selected Q matrix ordered by metadata. Label boundary optima and state the haploid, linked-marker limitation on the figure.

## Failure behavior

The figure stage refuses missing upstream states, input/state checksum disagreements, sample-to-metadata mismatches, malformed trees, malformed Q matrices, unknown taxa, non-finite coordinates where finite values are required, and pre-existing unvalidated outputs. It never silently drops samples or substitutes archived inputs.

## Testing and acceptance

Unit tests cover the fixed palette, signed-zero-centered FST normalization, deterministic distance-aware network positions, pie composition, unrooted tree layout, and required output formats. Workflow tests prove that `figures` occurs after ADMIXTURE and before reports, and that reports depend on the figure state. A miniature renderer test produces nonempty PNG/PDF/SVG files and a checksummed manifest. Layout tests continue to reject archive dependencies and absolute workstation paths.

The final acceptance gate requires a complete figure state, all declared figure files, exact checksums, no archive path in figure provenance, and a figure manifest included in the canonical deliverables manifest.
