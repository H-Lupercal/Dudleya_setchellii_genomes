# Dudleya supplementary analysis

This workspace implements approved decision plan v2.5 without modifying the canonical publication run. The immutable base is `canonical_publication/` run `publication-20260817`; all supplementary outputs belong to run `supplement-20260824` and remain under this directory.

## Supported command

```bash
supplementary_analysis/run_pipeline.sh \
  --config supplementary_analysis/config/supplementary_config.toml \
  --run-id supplement-20260824 \
  [--resume]
```

The entrypoint validates the canonical filesystem before every stage and again after acceptance. A resume is allowed only when the configuration, source code, imported canonical pure modules, canonical snapshot, upstream fingerprints, and saved outputs still match.

## Dependency order

```text
canonical_guard → metadata → identity → sensitivity → claims → inheritance
→ phase1_gate → likelihood_mapping → comparative_analyses → figures
→ reports → acceptance → canonical_guard_final
```

No preprocessing or mapping command is part of this workflow. Permissive and strict calls, mitochondrial mask variants, and all descendants use the immutable canonical BAMs. Raw reads are read only for provider-MD5 revalidation and Mash identity sketches.

## Interpretation boundaries

- DUSE is retained in sample-level PCA, trees, QC, and genotype displays but excluded from population inference while its label is unresolved.
- Mash similarity is a screen and never confirms sample identity by itself. A negative screen does not prove biological independence.
- Index hopping is untestable without actual index sequences and demultiplexing metrics.
- Chloroplast and mitochondrial trees and comparisons are unrooted.
- ADMIXTURE remains a demoted sensitivity display for linked haploid markers.
- Residual NUMT/NUPT ambiguity cannot be excluded because no nuclear decoy is available.
- Geography is `not_run:no_approved_coordinates`.

`CURRENT_RUN` is created or updated only after supplementary acceptance passes. Intermediate files, BAM-derived depth caches, IQ-TREE checkpoints, and raw sketches are ignored beneath `work/`; final tables, figures, reports, configurations, and provenance manifests are tracked.

