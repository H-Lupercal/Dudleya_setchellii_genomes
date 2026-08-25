# Dudleya supplementary analysis

This workspace implements approved decision plan v2.6 without modifying the canonical publication run. The immutable base is `canonical_publication/` run `publication-20260817`. The current supplementary run is `supplement-20260824-v26`; accepted v2.5 run `supplement-20260824` is preserved as superseded evidence.

The approved scope and scientific decision rules are recorded in [`supplementary_analysis_decision_plan.md`](supplementary_analysis_decision_plan.md).

## Supported command

```bash
supplementary_analysis/run_pipeline.sh \
  --config supplementary_analysis/config/supplementary_config.v2.6.toml \
  --run-id supplement-20260824-v26 \
  [--resume] [--from-stage STAGE]
```

Every executed stage validates the canonical filesystem. With `--from-stage --resume`, every preceding stage is actually executed in resume-validation mode before the requested stage, so changed fingerprints or outputs stop the run. The initial and final guards also verify the accepted v2.5 run-specific artifacts.

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
- PC–technical-covariate correlations cannot distinguish genuine divergence from reference-mapping bias.
- The mitochondrial mask-restricted likelihood map is diagnostic only and cannot trigger NeighborNet.

`status: PASS` means **workflow/provenance acceptance: PASS**; scientific claims retain separate PASS, PASS_WITH_CAVEAT, or FAIL statuses in the claim matrix. `CURRENT_RUN` is updated only after workflow/provenance acceptance passes. Intermediate files, BAM-derived depth caches, IQ-TREE checkpoints, and raw sketches are ignored beneath `work/`; final tables, figures, reports, configurations, and provenance manifests are tracked.
