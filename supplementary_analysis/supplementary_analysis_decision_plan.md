# Dudleya Organelle Population Genomics — Supplementary Analysis Decision Plan

**Version 2.6 · 2026-08-24 — FINAL POLISHING REVISION**

This revision supersedes v2.5 without declaring its accepted results faulty. The complete v2.5
plan is preserved at `decision_plans/supplementary_analysis_decision_plan.v2.5.md`, and its
accepted run remains `supplement-20260824` at Git commit
`dfbf23dacea8edd220cab88d090692e2cf7a5099`.

The current run is `supplement-20260824-v26`. It retains the six approved figure families and
all v2.5 scope exclusions. It adds only:

1. A formal claim for technical/reference-covariate associations, with wording that correlation
   cannot distinguish genuine divergence from reference-mapping bias.
2. A fully called-site PCA sensitivity using 1,111 chloroplast and 31 mitochondrial MAC≥2 SNPs,
   9,999 Procrustes permutations, seeds 424318–424319, and technical-covariate seeds 424320–424337.
3. A table of the ten largest π and FST sensitivity changes per nonzero scenario comparison,
   without changing the preregistered global PASS thresholds.
4. A diagnostic-only mitochondrial likelihood map restricted to the exact 43,182-base canonical
   high-confidence mask, using 100,000 quartets, `TPM3u+F+I`, and seed 314159.
5. Explicit separation of workflow/provenance acceptance from scientific claim statuses.
6. Real upstream fingerprint and output validation for `--from-stage --resume`.

The mitochondrial mask-restricted likelihood map does not replace the primary likelihood map and
cannot trigger NeighborNet. No mapping, NeighborNet, geography, STRUCTURE, UMAP/t-SNE, SFS,
mutation-spectrum, additional ADMIXTURE, or seventh figure-family analysis is authorized.

`canonical_publication/` and all v2.5 run-specific artifacts are immutable. `CURRENT_RUN` may be
updated to `supplement-20260824-v26` only after complete workflow/provenance acceptance passes.
