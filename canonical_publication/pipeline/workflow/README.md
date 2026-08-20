# Canonical dependency graph

This directory documents orchestration; it is not a second entrypoint. Run the
workflow only through `canonical_publication/run_pipeline.sh`.

The enforced dependency graph is:

```text
immutable sources ──┬──> references, annotations, static masks ──┐
                    └──> sample metadata ────────────────────────┼──> preprocessing + mapping
                                                               └──> mapping provenance gate
references + mapping provenance ───────────────────────────────────> QC eligibility + read-backed mtDNA mask
QC ──> all-site calls + filtered SNP layers ──┬──> callable consensus ──┬──> haplotypes
                                              │                        ├──> population genetics ──> independent cross-check
                                              │                        └──> primary organelle trees
                                              │                              ├──> supplementary partitioned concatenation
                                              │                              ├──> organelle-conflict analysis
                                              │                              └──> fixed-seed reproducibility
                                              ├──> organelle-specific PCA
                                              └──> supplementary ADMIXTURE
all accepted stages ──────────────────────────────────────────────────> figures + reports + invalidation + manifests
```

Every stage state stores its input and upstream fingerprints. `--resume`
validates those fingerprints and all declared output checksums; it raises a
stale-output error instead of silently reusing changed work.
