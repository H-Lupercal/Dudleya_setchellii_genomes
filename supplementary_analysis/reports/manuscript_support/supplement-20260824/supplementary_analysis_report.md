# Supplementary analysis report — supplement-20260824

## What we studied

We tested whether the canonical chloroplast and mitochondrial conclusions are sensitive to
approved filtering and mitochondrial-mask choices. We quantified phylogenetic information,
compared supported unrooted organelle topologies, examined technical covariates, and
standardized marker and sample counts.

## Why it matters

Organelle genomes are linked haploid lineages. Strong-looking clusters or differentiation can
reflect filtering, missingness, reference concordance, or limited phylogenetic information.
These analyses expose those alternatives without treating organelles as independent nuclear loci.

## Data and scope

The immutable base is `publication-20260817` (cp 276; mt 271; shared 271). Existing filtered
BAMs were reused. No preprocessing or remapping was performed. DUSE samples remain in
sample-level displays but are excluded from population inference because the population label
is unresolved. Geography was `not_run:no_approved_coordinates`.

## Robustness outcomes

- permissive chloroplast pi: PASS
- permissive chloroplast fst: PASS
- permissive chloroplast pca: PASS
- permissive mitochondria pi: PASS
- permissive mitochondria fst: PASS
- permissive mitochondria pca: PASS
- strict chloroplast pi: PASS
- strict chloroplast fst: PASS
- strict chloroplast pca: PASS
- strict mitochondria pi: PASS_WITH_CAVEAT
- strict mitochondria fst: PASS
- strict mitochondria pca: PASS
- mtmask70 mitochondria pi: PASS
- mtmask70 mitochondria fst: PASS
- mtmask70 mitochondria pca: PASS
- mtmask90 mitochondria pi: PASS
- mtmask90 mitochondria fst: PASS
- mtmask90 mitochondria pca: PASS

Scientific caveats or failures remain visible in the figures and require interpretation changes
recorded in `claim_analysis_decisions.tsv`; they are not provenance failures.

## Phylogenetic information

- chloroplast: resolved=98.11%, partly resolved=0.59%, unresolved=1.30% (TREE_LIKE_NO_NETWORK); composition failures=1/276, >50% gaps/ambiguity=0.
- mitochondria: resolved=82.01%, partly resolved=3.53%, unresolved=14.46% (TREE_LIKE_NO_NETWORK); composition failures=269/271, >50% gaps/ambiguity=271.

The cp–mt comparison is unrooted. After contracting branches lacking joint SH-aLRT≥80 and
UFBoot≥95 support, normalized RF was 64/96 =
0.6667 on 229 mitochondrial unique-sequence representatives.
The tanglegram restores/displays all 271 shared samples, including 42 annotated zero-length
identical-tip memberships.

Across eligible common population pairs, cp–mt FST rank agreement was
rho=0.7956 using
557/561 finite pairs;
4 nonfinite pairs remain in the table but were not
used in the correlation.

## Marker-count and sample-size controls

Observed mitochondrial multi-population haplotype sharing was
11; this equals the median of the
1,000 chloroplast 146-site draws (95% interval
8–
15). Thus this sharing count is
marker-count-sensitive. At common n=4, CY_SIE and CY_CAS retained the two highest median
chloroplast pi values; their medians exceeded
100.0% and
96.7% of draws from other
populations, respectively. This is a descriptive resampling comparison, not an independent-locus
hypothesis test.

## Interpretation

The organelles describe lineage history, not a nuclear-genome admixture history. ADMIXTURE
remains a demoted sensitivity visualization because its linked haploid markers violate the usual
independent-diploid interpretation. Organelle inheritance mode was not established here.

## Limitations

- No nuclear decoy was available, so residual NUMT/NUPT ambiguity cannot be excluded.
- Raw-read sketches are screens: sketch similarity alone is suspected, not confirmed identity,
  while a negative result does not prove biological independence.
- Index hopping is untestable without index sequences, sample sheets, and demultiplexing metrics.
- Mitochondrial inference receives extra restraint because the primary alignment contains only
  146 SNPs, 269/271
  sequences failed the composition test, and weak-split topology was not fully reproducible.
- Marker-count resampling does not control mutation rate, mask, missingness, or organelle biology.
