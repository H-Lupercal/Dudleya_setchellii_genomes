# Canonical publication analysis

This tree contains the repaired Dudleya organelle workflow and only outputs
whose complete dependency chain was generated or fingerprint-validated here.

Create the pinned main environment at the runner's local tool prefix (the
prefix is ignored by Git), and create the independent validator separately:

```bash
micromamba create -y -p .tools/bioconda-env \
  -f canonical_publication/environment.yml
micromamba create -y -p .tools/scikit-allel-validation \
  -f canonical_publication/validation_environment.yml
```

Run from the repository root:

```bash
canonical_publication/run_pipeline.sh \
  --config canonical_publication/config/publication_config.toml \
  --run-id YYYYMMDD-label
```

The chloroplast and mitochondrial analyses use independent QC-derived sample
sets. ADMIXTURE and the partitioned concatenated tree are supplementary.
Because no nuclear assembly is available, NUMT/NUPT ambiguity is reduced by
MAPQ, base-quality, depth, and duplicate filters but cannot be eliminated.

The independent population-genetics acceptance check uses the separately
pinned `validation_environment.yml`; the main environment is pinned in
`environment.yml`. The runner defaults to the validator prefix shown above;
`--validation-python` is available only when a different separately pinned
validator is required. Add `--resume` only to validate and reuse outputs whose
complete fingerprints still match.
Stage slicing with `--from-stage` is diagnostic/operational only: it requires
`--resume`, and the runner revalidates every preceding stage before executing
the requested downstream stage.

Variant products are deliberately layered: `high_confidence_variant_sites`
retains fixed-alternate sites needed for correct consensus sequences,
`primary` contains segregating sites including singletons, and
`mac2_ordination` is restricted to MAC≥2 for PCA and supplementary ADMIXTURE.
Raw mpileup likelihood BCFs are also retained as ignored work products because
bcftools does not emit GQ at invariant `ALT=.` calls. Callable consensus
generation derives haploid invariant GQ from the homozygous reference and
symbolic-nonreference PLs, while accepted biallelic SNPs use the masked haploid
bcftools GT/GQ calls. Both paths enforce DP≥5 and GQ≥20; uncertain positions
remain `N`.

All 278 complete read pairs in this immutable deposit contain one balanced
lane. The streaming mapper fails explicitly on a future multi-lane manifest
instead of silently concatenating lanes without validated lane-aware
provenance.

Canonical scientific stages refuse inputs resolving beneath
`archive_noncanonical/`. The final audit-report stage reads the preserved
manifest from `canonical_publication/provenance/archive/` to account for every
retired artifact without requiring the archived snapshot in a canonical
checkout. Tracked snapshot contents remain available on
`archive/noncanonical-2026-08-17`; archived scientific outputs are never
analysis inputs.
