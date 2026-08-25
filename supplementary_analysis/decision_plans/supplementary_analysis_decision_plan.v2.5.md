# Dudleya Organelle Population Genomics — Supplementary Analysis Decision Plan

**Version 2.5 · 2026-08-23 — FINAL, APPROVED**
**Supersedes:** v1.0 master catalogue; v2.0 first reduction; v2.1 scientific corrections; v2.2 execution definitions; v2.3 approval corrections; v2.4 π/reproducibility corrections.
**This version** fixes the one remaining defect — the v2.4 NeighborNet trigger was logically impossible (two "incompatible" splits within one inferred tree cannot exist; the referenced conflict table is cross-organelle). The trigger now uses within-organelle **bootstrap split frequencies**. Scope is final: **4 verification checks + 6 new figure families + presentation replacements; geography and NeighborNet conditional.** No analysis family is added or removed. Every repo-fact below was re-checked read-only this session.

The canonical run being verified/extended is `publication-20260817` (acceptance PASS). Nothing here revives the pre-remediation archive logic.

> **Path convention:** all paths are from the repository root and begin with `canonical_publication/`.
> **Documentation location:** `supplementary_analysis/supplementary_analysis_decision_plan.md`. Planning artifact only; no code or canonical output modified.

---

## 0. Review-response changelog

### v2.4 → v2.5 (sixth review — 1 logical correction; approved, +1 implementation note)

| # | v2.4 issue | Fix in v2.5 | § |
|---|---|---|---|
| 1 | NeighborNet trigger required "≥1 pair of incompatible splits each SH-aLRT≥80 & UFBoot≥95 **within** that organelle" — **impossible**: all splits in one inferred tree are mutually compatible, and `strongly_supported_organelle_conflicts.tsv` is a **cross-organelle** (cp-split vs mt-split) table, not within-organelle. | Trigger now measures within-organelle conflict from the **1,000 UFBoot bootstrap replicate splits** (`{organelle}.primary.splits.nex`): ≥1 pair of **mutually incompatible splits each recovered in ≥20% of replicates**. Branch SH-aLRT/UFBoot labels are not used for this; the cross-organelle table is correctly assigned to Figure 4. | 5, 6 |
| 2 | (approval implementation note) mt `.splits.nex` has **229** unique taxa, not 271 — IQ-TREE collapsed 42 identical mt sequences (cp uncollapsed at 276). | Captured: the ≥20% criterion runs over the 229-taxon mt split space; mt tree/network presentation **restores/annotates identical samples** via `mitochondria.primary.uniqueseq.phy` / log NOTEs. Not a plan defect. | 5, 6 |

### v2.3 → v2.4 (fifth review — 1 required correction + 3 ambiguities + 1 clarification)

| # | v2.3 issue | Fix in v2.4 | § |
|---|---|---|---|
| 1 | Fig 5 π resampling used "confidently-assigned samples only" | **Corrected — π uses all QC-eligible samples**, recomputing the jointly-callable denominator per draw (haplotype ambiguity reflects missing positions, not absence of callable π; verified: ABAB_MAD is fully haplotype-ambiguous yet π=0.000214). Confidently-assigned eligibility applies only to haplotype richness/sharing. Audit-table diversity row split accordingly. | 6 |
| 2 | likelihood-mapping input "or a documented derivative" | **Committed** to the canonical callable alignment. | 5 |
| 3 | NeighborNet trigger "high side fraction" was subjective | **Numerical + evidence-based**: side fraction > 20% **AND** ≥1 pair of incompatible splits each SH-aLRT≥80 & UFBoot≥95 (supported conflict, the `strongly_supported_organelle_conflicts.tsv` criterion). A high side fraction alone does not trigger it. | 5, 6 |
| 4 | cp downsampling "many (e.g. 1,000)" non-deterministic | **Exactly 1,000 replicates, fixed seed 424200**; π resampling fixed seed 424201. | 6 |
| 5 | coordinates added to byte-frozen canonical `populations.tsv` | Future coordinates **versioned under the supplementary run**, never the frozen canonical file. | 3.5, 6 |

### v2.2 → v2.3 (fourth review approval — 3 required corrections + 2 clarifications)

| # | v2.2 issue | Fix in v2.3 | § |
|---|---|---|---|
| 1 | NeighborNet triggered on *central/unresolved* signal | **Trigger logic corrected**: central quartets = insufficient information (a network cannot recover absent signal → no auto-NeighborNet); NeighborNet only for high *partly-resolved/conflicting* signal or contradictory well-supported patterns. | 5, 6 |
| 2 | likelihood-mapping "fixed seed" not a value | **Predeclared** cp `-seed 271828`, mt `-seed 314159` (the canonical `cp_seed`/`mt_seed`); model citations moved to the **primary** `.iqtree` files. | 5 |
| 3 | claim-matrix & inheritance safeguards unassigned | Added as **Phase 1 gating documentation deliverables** with named artifacts (`claim_analysis_decisions.tsv`, `organelle_inheritance_evidence.md`). | 6 |
| 4 | "two full recalls" ambiguous | Clarified: **re-genotype from fingerprint-validated existing BAMs/pileups + downstream — no remapping of raw reads.** | 6 |
| 5 | raw-read sketch could imply proof | Sketches are duplicate-library **screens** with predeclared thresholds/known controls; a negative sketch alone does **not** establish field-sample independence. | 3.2 |

### v2.1 → v2.2 (third-round follow-up — execution definitions + 2 method fixes)

| # | v2.1 issue | Fix in v2.2 | § |
|---|---|---|---|
| 1 | DUSE correction targeted `populations.corrected` | **Authoritative correction is sample-level** (`samples.corrected-*.tsv`); population table is *derived* from it. Correction manifest schema defined. | 2.5 |
| 2 | `S###`+lane called "exactly what an index-hopping check needs" | `S244` is a demux sample number, **not the index sequence**. Index-hopping needs index seqs/sample sheet/demux metrics → else **untestable**. Duplicate-library detection uses **raw-read sketches**, not organelle genotype identity. | 3.2 |
| 3 | Sensitivity used rank correlation only | Added **magnitude criteria** (proportional/absolute change); ρ=1 can hide a 10× shift. | 6 |
| 4 | Sample-count row required 276/271 under all scenarios | **Split**: canonical baseline exact 276/271; sensitivity counts may change per eligibility rule and are reported, not failed. | 6 |
| 5 | `unresolved` identity outcome had no next step | **Four outcomes** (verified / confirmed / suspected / unresolved) each with an action. | 3.2, 6 |
| 6 | "one mt-mask-only run"; 3 full re-calls implied | Corrected workload: **2 new full recalls + 2 mt-mask reruns**; canonical reused after fingerprint check (80% mask and canonical scenario already exist). | 6 |
| 7 | Likelihood mapping had no reproducibility/escalation spec | Added quartet count, seed, model, input, and a **numeric "substantial unresolved" threshold** driving the NeighborNet decision. | 5 |
| 8 | Figure-5 resolution matching underdefined | **Repeated seeded downsampling** distribution; framed as marker-count sensitivity only; **DUSE excluded from resampling**; common n declared. | 6 |
| 9 | "H23 is a DUSE artifact" | Softened: the **haplotype is real**; its **cross-population classification** is DUSE-label-dependent. | 2.1 |
| 10 | RF/PCA/TUL2 execution gaps | RF label-parsing + multifurcation support; PCA permutation count for p<0.001; TUL2 display-name resolution added to Phase 1. | 6 |

### v2.0 → v2.1 (third review — scientific corrections)

| # | v2.0 claim | Correction | Verified against |
|---|---|---|---|
| 1 | DUSE+CY_CAS+CY_SIE are input defects | **Only DUSE** is a label defect; CY_CAS = species-field inference; CY_SIE = statistical outlier, not a metadata entry | `...source_metadata_ambiguities.tsv` |
| 2 | "DUSE in 34 of 595 pairs" as evidence | Removed — every pop is in exactly 34 pairs | arithmetic |
| 3 | 0/50 vs 6/45 = "corrected result" | Relabelled DUSE-excluded sensitivity | 2.1 |
| 4 | sCF "undefined" for one molecule | sCF is **site-based, computable**; dropped for redundancy | IQ-TREE docs |
| 5 | reproducibility table "better than likelihood mapping" | Contradiction removed; measures weak-split instability, not information content | `...fixed_seed_reproducibility.tsv` |
| 6 | identity audit "rules out" | Audit **flags**; adds provider MD5; reports `unresolved` | `...provider_md5_validation.tsv` |
| 7 | "five of six figures need no compute" | False; corrected accounting | own table |
| 8 | multiple-testing on 595 FST | FST are effect-size estimates with CIs; correction only where p-values exist | `...pairwise_hudson_fst.tsv` cols 10–11 |

---

## 1. Non-negotiable baseline and interpretation rules

| Canonical element | State |
|---|---|
| Canonical run | `publication-20260817`, acceptance PASS |
| Chloroplast / Mitochondrial cohort | 276 / 271 QC-eligible samples |
| Shared cohort | 271 (cp∩mt) — direct cp-vs-mt comparison and concatenation only |
| Populations | 35 population codes |
| Primary segregating SNPs | cp 2,261 · mt 146 |
| High-confidence sites (incl. fixed-alt) | cp 2,273 · mt 157 |
| Variant layers | high-confidence (keeps fixed-alt) · primary MAC≥1 incl. singletons · ordination/ADMIXTURE MAC≥2 |
| QC thresholds | haploid; DP≥5; GQ≥20; site QUAL≥30; ≤20% missingness; eligibility ≥80% breadth at DP≥5 over the organelle unique-mappability denominator |
| mt high-confidence mask | read-backed; built at ≥80% eligible-sample support |
| Phylogeny | separate cp & mt **unrooted** IQ-TREE; 1,000 SH-aLRT + 1,000 UFBoot with BNNI |
| ADMIXTURE | K=1–12, 10 seeded replicates/K; min mean CV at K=12 (upper boundary tested) |

**Interpretive guardrails:** (1) organelle, not nuclear — cp/mt SNPs are effectively haploid and highly linked; (2) ADMIXTURE/STRUCTURE bars are organelle haplotype-clustering/sensitivity, never ancestry/hybridization/introgression; (3) mt inference uses the restricted read-backed mask, not the 243,359 bp candidate; (4) NUMT/NUPT flagged, never eliminated; (5) trees stay **unrooted** — external refs (*D. farinosa* NC_085682.1; *G. paraguayense* PV256627.1) are not automatically valid outgroups; (6) say "identical organelle haplotypes," never "clonality"; (7) mt gets extra restraint (146 SNPs).

---

## 2. Evidence (re-checked read-only)

### 2.1 BLOCKER — DUSE is the only confirmed population-label defect

`canonical_publication/metadata/qc/publication-20260817/source_metadata_ambiguities.tsv` has exactly four entries: `CY_BOU`, `CY_CAS` (blank **species** inferred; population labels intact), `DUSE` ("population code declared only in source column header"), `TUL2` (conflicting source population labels). **CY_SIE is absent.**

- **DUSE — population-label blocker.** Resolve first (§2.5).
- **CY_CAS — taxonomy check only.** Verify inferred *D. cymosa*; label is fine.
- **CY_SIE — biological/statistical outlier to verify**, not a presumed metadata error.

High π + low mean FST can indicate heterogeneous lineages but does not prove a label is wrong. Profile (`chloroplast.population_summary.tsv` + `chloroplast.pairwise_hudson_fst.tsv`):

| pop | cp π | mean cp FST | status |
|---|---|---|---|
| CY_SIE | 0.00291 | 0.664 | outlier to verify (not in ambiguity table) |
| CY_CAS | 0.00211 | 0.546 | taxonomy to verify; label fine |
| DUSE | 0.00205 | 0.449 | **label blocker** |
| all others | ≤0.00043 | ≥0.786 | |

**Haplotype-sharing contrast — DUSE-excluded sensitivity, not "corrected":**

| | multi-pop haplotypes | DUSE excluded |
|---|---|---|
| cp | 1 of 50 | 0 of 50 |
| mt | 11 of 45 | 6 of 45 |

`H23 = DUSE:1, TULP:2` **is a real haplotype**; its *classification as cross-population* is DUSE-label-dependent, and 5 of 11 mt shared haplotypes involve DUSE. Keep descriptive and DUSE-excluded until DUSE is resolved and the comparison is resolution-matched (§6-Fig5).

### 2.2 STRIKE — independent π and FST already done
`...trusted_scikit_allel_crosscheck.tsv`: 4,830 comparisons, all match, max deviation 5×10⁻¹⁰. Pairwise FST already carries `bootstrap_ci_2.5`/`bootstrap_ci_97.5`. Cite; do not rebuild.

### 2.3 Tree reproducibility — scope of the claim
`...fixed_seed_reproducibility.tsv`: cp RF=0 (90 strong splits); **mt RF=64** of 268 internal (17 strong splits, all reproduced); concatenated RF=0. This establishes **weak-split topology instability** in mt, **not** alignment information content (that is likelihood mapping, §5). Operational consequence: any cp-vs-mt RF must contract weak branches first (§6-Fig4).

### 2.4 Patterns worth showing (verified)
π=0 populations: cp 7 (KIRT n=10, QUIN n=10, QUI1 n=9, …), mt 13. FST saturation: 74% of cp pairs >0.9, median 0.964, one negative overall. Missingness: 180/276 cp and 72/271 mt `AMBIGUOUS`. cp SNP density per 5 kb: 60–130 background, 163 & 173 at 140–150 kb, zero blocks at 75–95 & 115–140 kb. cp-vs-mt FST concordance r=0.93, ρ=0.82 (591 finite pairs). PCA PC3: cp 12.3%, mt 10.6%.

### 2.5 DUSE correction — sample-level provenance & invalidation policy

Reassigning DUSE means changing each affected sample's `popcode`. `canonical_publication/metadata/populations/populations.tsv` is only popcode→species→name; the sample→popcode map lives in `canonical_publication/metadata/samples/samples.tsv`. Therefore:

1. **Authoritative correction is sample-level:** write `canonical_publication/metadata/samples/samples.corrected-YYYYMMDD.tsv`. **Derive** a regenerated population table from it; do not hand-edit the population table.
2. **Correction manifest** (one row per changed sample): `sample_id`, `old_popcode`, `new_popcode_or_EXCLUDED`, `evidence_source`, `decision_author`, `decision_date`, `confidence_or_unresolved`.
3. Allocate a **new supplementary run ID** (`supplement-YYYYMMDD`); leave `publication-20260817` **byte-unchanged**.
4. **Regenerate every population-dependent output** under the new run ID; never edit canonical files in place.
5. **Manifest** links each replacement output to the corrected-metadata fingerprint.
6. **Fallback if source records cannot resolve DUSE:** exclude DUSE from *population-level* inference while **retaining its samples in sample-level analyses** (PCA points, tree tips, per-sample QC). Record as an explicit cited decision, not a silent drop.

---

## 3. Missing safeguards

### 3.1 Claim→analysis matrix + stopping rules
Every analysis names the claim it supports and the result that would change the interpretation; if it cannot change a conclusion, delete it. It also prunes Figure-5 components with no attached claim. Predeclare the negative-result rule (e.g. "threshold-dependent structure → report ranges, not clusters").

### 3.2 Sample-identity & independence audit — flags, does not rule out

Evidence, strongest first:
1. **Exact duplicated files** — provider MD5 in `...publication-20260817.provider_md5_validation.tsv`.
2. **Structured IDs** — specimen/extraction (`Du-###`), plate/well (`LP_###`). **Index hopping:** the filename `S###` is a *demultiplexed sample number, not an index sequence*; a defensible test needs index sequences / the sample sheet / demultiplexing metrics / unexpected index combinations. **If those are unavailable, mark index hopping `untestable`** — do not infer it from `S###`+lane.
3. **Duplicate libraries:** use **raw-read (or nuclear-read) genome sketches** (e.g. MinHash/`mash`) as a **screen** with predeclared similarity thresholds and at least one known-control pair (a same-library and a different-library pair, if identifiable). A negative sketch result alone does **not** establish field-sample independence. **Do not** infer duplicate libraries from identical *organelle* genotypes — unrelated samples can share an organelle haplotype.
4. **Contamination:** haploid mixed-allele fractions and depth patterns where raw pileups permit.

**Four outcomes, each with an action:**
- **Verified independent** → proceed.
- **Confirmed duplicate/mislabel** → correct or exclude, and invalidate the downstream results built on it (§2.5 machinery).
- **Suspected** → run the affected analysis with-and-without, as a sensitivity check.
- **Unresolved** (records absent) → retain only with an explicit stated limitation; **do not claim biological replication or population clonality.**

Outcome language throughout: "identical organelle haplotypes," never "clonality."

### 3.3 Uncertainty & multiplicity policy
Descriptive pairwise FST uses the existing bootstrap CIs and effect sizes — **not** p-value multiplicity correction. Apply correction **only where p-values are generated** (AMOVA permutations, PCA–QC correlation tests, Procrustes/protest). Predeclare CI method, permutation counts, minimum usable population size, and the correction family per test. mt claims get extra restraint.

### 3.4 Organelle-inheritance justification
Verify or cite maternal inheritance of cp *and* mt in *Dudleya*/Crassulaceae; if unsupported, soften inheritance language throughout.

### 3.5 Sensitive-location policy
Any coordinates added later need provenance and permissions; public outputs generalize or withhold exact localities. Coordinates are **versioned under the supplementary run** (e.g. `canonical_publication/metadata/populations/populations.coords-<run-id>.tsv`), **never** added to the byte-frozen canonical `populations.tsv`. Precondition on §7 geographic work.

---

## 4. Merged cut list

**Keep:** scoped likelihood-mapping panel (§5); tanglegram + one support-filtered RF + FST scatter; one integrated missingness/confounder figure; one targeted resampling; genome-coordinate track; PCA scree + PC3.

**Cut:** sCF (computable but **redundant** with UFBoot/SH-aLRT, and linked sites aren't independent replicates — **not** "undefined"); Procrustes + quartet distance; mutation spectrum / substitution classes / SFS; per-gene burden; NUMT/NUPT as its own figure (→ flag columns in Fig 3); the full filtering grid; STRUCTURE; extended-K and native-haploid ADMIXTURE reruns; DAPC/UMAP/t-SNE; cloudograms/271-tip bootstrap views; rooted trees; species-level AMOVA unless the hierarchy/sampling make it interpretable (within-*D. setchellii* is more defensible — 29/35 populations are *setchellii*). **Ti:Tv** removed as a deliverable (compute only inside variant-call validation if a specific concern arises).

**Conditional:** NeighborNet (§5); geography (§3.5, §6).

---

## 5. Likelihood-mapping decision (with reproducibility + escalation)

- **sCF: cut** (§4).
- **Likelihood mapping: keep, one panel per organelle** — it answers whether the quartet signal is tree-like at all, which the reproducibility table does not. Predeclared, reproducible spec:
  - **Input (committed):** the canonical callable alignment per organelle — `canonical_publication/results/alignments/publication-20260817/{organelle}.callable_alignment.fa`. No derivative.
  - **Quartets:** 100,000 random quartets (`-lmap 100000`).
  - **Seed (predeclared):** cp `-seed 271828`, mt `-seed 314159` — the canonical `cp_seed`/`mt_seed` from `canonical_publication/config/publication_config.toml`, so the mapping matches the primary-tree runs.
  - **Model:** reuse the canonical ModelFinder model per organelle — cp `TVM+F+I+R4` (`canonical_publication/results/trees/publication-20260817/chloroplast.primary.iqtree:44`); mt `TPM3u+F+I` (`.../mitochondria.primary.iqtree:44`).
  - **Report:** the 7-region percentages; "resolved" = 3 corners, "partly/conflicting" = 3 sides, "unresolved/uninformative" = center.
  - **Decision (central signal is *absence* of information, which a network cannot recover; a high side fraction is partial resolution, not proof of conflict):**
    - **High central/unresolved fraction** (center **> 15%**) → report **insufficient phylogenetic resolution** for that organelle; **do not** run NeighborNet.
    - **NeighborNet runs only on demonstrated within-organelle conflict** — both conditions: (a) side/partly-resolved fraction **> 20%**, **AND** (b) the organelle's own **bootstrap replicate splits** — from the 1,000 UFBoot trees, as recorded in `canonical_publication/results/trees/publication-20260817/{organelle}.primary.splits.nex` — contain **≥ 1 pair of mutually incompatible splits each recovered in ≥ 20% of replicates**. This measures competing, well-supported alternative resolutions, which is exactly what a split network displays. It deliberately does **not** use branch SH-aLRT/UFBoot labels — those attach to a single tree whose splits are all mutually compatible — and is unrelated to `strongly_supported_organelle_conflicts.tsv`, which is a **cross-organelle** (cp-split vs mt-split) table used in Figure 4, not a within-organelle one.
    - **Clearly tree-like** (resolved corners dominate, no competing supported splits) → no NeighborNet.
  - **Implementation note (taxon collapse):** cp `.splits.nex` spans all **276** taxa, but mt `.splits.nex`/UFBoot span **229 unique taxa** — IQ-TREE collapsed 42 identical mt sequences (log NOTES "X is identical to Y but kept…"; mapping also in `mitochondria.primary.uniqueseq.phy`). The ≥20%-replicate criterion operates over the 229-taxon mt split space; that is correct for detecting conflict, but any mt **network or tree presentation must restore/annotate the 42 identical samples** so none is silently dropped (see Figure 4).

---

## 6. Reduced plan — 4 checks + 6 figures (hard cap)

### Phase 1 — Verification (gates Phase 2)

1. **Resolve DUSE** (§2.1, §2.5) — BLOCKER, first. Separately **verify CY_CAS taxonomy** and **verify the CY_SIE outlier** (neither is a presumed label defect). **Resolve the TUL2 display name** (`"Dud Tulare Hill-2"` vs `TH Clone site 2`) before any figure labels it.
2. **Sample-identity & independence audit** (§3.2) with the four-outcome rule, including the π=0 duplicate check via raw-read sketches.
3. **Re-verify provenance + acceptance** — checksums, `canonical_publication/provenance/runs/publication-20260817/ACCEPTANCE.json`, internal consistency.
4. **Filtering sensitivity.** Scenarios and **actual workload**:
   - **canonical** — reuse existing outputs after fingerprint validation (**no recall**);
   - **permissive** (DP3, GQ15, 30% miss, breadth 70) — **new full recall**;
   - **strict** (DP10, GQ30, 10% miss, breadth 90) — **new full recall**;
   - **mt-mask-support** 70% and 90% — **two new mt-mask reruns** (mt only, downstream of the mask; 80% is canonical and reused).
   - **New compute total: 2 full recalls + 2 mt-mask reruns** — where "recall" means **re-genotyping from the fingerprint-validated existing BAMs/pileups and re-running the affected downstream stages (eligibility, variants, popgen, PCA, trees), NOT remapping raw reads.** Mapping/preprocessing is reused after fingerprint validation; there is no multi-billion-read rerun.
   - **Held fixed** (limitation, not varied): QUAL, MAPQ, base-quality. Bundled scenarios show *whether* results move, not *which* threshold moved them — an accepted trade.
   - **Metrics per scenario, across populations/pairs present and adequately sampled in both** (see thresholds): samples retained; SNPs retained; **π** — Spearman ρ **and** median/max proportional change vs canonical; **FST** — Spearman ρ **and** median/max absolute change; **PCA** — Procrustes/`protest` correlation (≥9,999 permutations so p<0.001 is attainable). **Predeclare** which numeric change would alter a headline claim.
5. **Claim→analysis decisions** (§3.1) — documentation gate, not an analysis. Produce `claim_analysis_decisions.tsv` (columns: `claim`, `analysis`, `metric`, `pass_caveat_fail`, `interpretation_change`). Any figure whose row carries no manuscript claim is dropped before Phase 2.
6. **Organelle-inheritance evidence** (§3.4) — documentation gate. Produce `organelle_inheritance_evidence.md`: the *Dudleya*/Crassulaceae inheritance evidence and the **exact manuscript language it permits**; if evidence is absent, inheritance wording is softened accordingly. (Both artifacts under the supplementary run's `reports/manuscript_support/<run-id>/`.)

**Struck as already-done:** independent π, independent Hudson FST (§2.2).

**Audit table — predeclared tolerances:**

| Component | Expected | PASS | PASS-with-caveat | FAIL |
|---|---|---|---|---|
| Provenance / checksums | exact | 100% exact | — | any mismatch |
| **Canonical baseline** sample counts | cp 276 / mt 271 | exact | — | any change |
| **Sensitivity** sample counts | vary by eligibility rule | reported + explained by the rule | — | change not explained by the rule |
| π under sensitivity | stable | ρ ≥ 0.95 **and** median proportional Δ ≤ 10% | 0.90 ≤ ρ < 0.95 or 10–25% Δ | ρ < 0.90 or median Δ > 25% |
| FST under sensitivity | stable | ρ ≥ 0.95 **and** median abs Δ ≤ 0.05 | 0.90 ≤ ρ < 0.95 or 0.05–0.10 | ρ < 0.90 or median abs Δ > 0.10 |
| PCA under sensitivity | stable | protest r ≥ 0.90 (p<0.001) | 0.80 ≤ r < 0.90 | r < 0.80 |
| Permissive vs strict | same conclusions | agree on all headline claims | agree on direction, differ in magnitude → report as range | disagree on a headline claim |
| Diversity claim — π/FST | adequately sampled | **eligible** n ≥ 5 | eligible n = 3–4 (flagged) | eligible n < 3 → no population-level claim |
| Diversity claim — haplotype richness/sharing | adequately sampled | **confidently-assigned** n ≥ 5 | assigned n = 3–4 (flagged) | assigned n < 3 → no haplotype claim |
| Identity/independence | resolved | verified independent | suspected → with/without sensitivity | confirmed duplicate/mislabel → correct/exclude+invalidate; unresolved → retain with limitation, no replication/clonality claim |

### Phase 2 — Supplementary figures (**6 named, hard maximum**)

| # | Figure | Answers | Computation |
|---|---|---|---|
| 1 | **Robustness** — 3 scenarios + mt-mask 70/90: π/FST ρ+magnitude, SNP & sample counts, PCA protest | Do conclusions depend on thresholds? | 2 full recalls + 2 mt-mask reruns + protest |
| 2 | **Phylogenetic information** — cp & mt likelihood-mapping ternary + region fractions (§5) | Is the (esp. mt) alignment tree-like? | IQ-TREE `-lmap 100000` (light) |
| 3 | **Genotype & technical-confounder** — sample×SNP + callability heatmaps + missingness/depth/reference-identity/NUMT-NUPT flag tests | What drives clusters; is missingness/reference bias faking structure? | heatmaps + correlation tests + flag columns |
| 4 | **cp–mt comparison** — tanglegram (271 shared) + support-filtered normalized RF + cp-vs-mt FST scatter | Do the histories agree, and where not? | prune + contract + RF + join |
| 5 | **Population diversity** — sampling-standardized diversity/haplotype sharing + π-vs-mean-FST scatter with sample sizes | Genuinely diverse/isolated vs undersampled? | targeted resampling |
| 6 | **Genome-coordinate track** — SNP density + callable/repeat/IR masks + projected annotations + depth, per callable kb | Where does callable variation sit? | window counts + normalization |

**Figure-4 RF rule:** prune both trees to shared samples; parse the IQ-TREE internal-node label as `SHaLRT/UFBoot`, contracting a branch unless **SH-aLRT ≥ 80 AND UFBoot ≥ 95** (the canonical "strong" definition); use an RF implementation that accepts **multifurcating** contracted trees (e.g. ete3/dendropy, unrooted); report **normalized RF with its denominator**; describe as **supported-topology compatibility**, not total evolutionary disagreement. **Taxon-collapse handling:** the mt tree/splits carry 229 unique taxa (42 identical mt sequences collapsed; cp is uncollapsed at 276). Before the tanglegram/RF, **restore the identical mt samples from the `mitochondria.primary.uniqueseq.phy` / log-NOTE mapping** so all shared samples are represented (identical samples attach to the same tip, annotated), and state whether RF is computed on the 229-unique or the restored tip set.

**Figure-5 rules:**
- (a) **Resolution matching by repeated seeded downsampling** — draw **exactly 1,000 subsets** of 146 cp sites (fixed **seed 424200**) on the shared 271 samples and compare observed mt sharing against the resulting cp *distribution*; describe as **marker-count sensitivity only**, not a control for mutation rate, mask, missingness, or genome biology.
- (b) **π resampling uses ALL QC-eligible samples**, recomputing the jointly-callable denominator for every draw — **not** confidently-assigned-haplotype samples. Haplotype ambiguity primarily reflects missing positions; ambiguous samples still contribute valid callable-site π (e.g. ABAB_MAD is fully haplotype-ambiguous yet has π = 0.000214). Excluding them — 180/276 cp — would bias diversity. **Exactly 1,000 draws, fixed seed 424201.** The diversity panel **excludes DUSE while its grouping is unresolved**; CY_SIE (eligible n=4) and CY_CAS (eligible n=5) bind the common size → **common n=4**; compare each outlier's π against same-n draws from other populations.
- (c) **Haplotype richness/sharing uses confidently-assigned haplotypes only**, with ambiguous counts reported separately. This eligibility rule is distinct from (b).

**Presentation fixes — folded into existing canonical figure replacements, NOT new supplement figures:** re-encode the FST heatmap (rank/quantile color, clustered/clade order, negatives visible, companion ECDF; numbers unchanged); PCA scree + PC1–PC3; population-collapsed annotated tree (strips: population, π, haplotype, callable fraction), full 271-tip tree kept archival; existing ADMIXTURE panel retained as one demoted sensitivity figure with the linked-marker caveat.

### Phase 3 — Public-facing write-up
**What we studied → Why it matters → What data → What we found → What it means → Limitations.** Each of π/FST/haplotype/PCA/bootstrap keeps a plain-English gloss. Limitations carry the organelle-vs-nuclear caveat, the ADMIXTURE demotion, the mt-resolution caveat, and any `unresolved` independence flags.

### Conditional (not counted in the 6)
- **Geographic / IBD** — only after coordinates are obtained, provenance/permission-checked, safe to disclose (§3.5), and tied to a pre-registered hypothesis. **No coordinates exist in the repo today**; verified `latitude`/`longitude`, **versioned under the supplementary run** (not the frozen canonical `populations.tsv`, per §3.5), unlock more than any plotting here.
- **NeighborNet** — only if Figure 2 meets **both** §5 conditions: likelihood-mapping side fraction > 20% **and** ≥1 pair of mutually incompatible splits each recovered in ≥20% of the organelle's 1,000 UFBoot replicates. A high *central/unresolved* fraction, or a high side fraction without competing supported bootstrap splits, does **not** trigger it.

---

## 7. Tool roles (reference)

FigTree/PearTree (inspect/style; iTOL interactive; **ggtree** preferred for scripted figures). Geneious/UGENE/Jalview/IGV = manual inspection, not a pipeline replacement. PopART (haplotype networks). SplitsTree (conditional; boxes descriptive, not hybridization proof). ADMIXTURE = one demoted panel. **Future nuclear data** (out of scope): mapping → GLs/called SNPs → filtering → relatedness removal → LD pruning → PCA/structure; NGSadmix/PCAngsd (low/uneven depth); conStruct/EEMS (spatial). Scale ≈ samples × retained SNPs, not reference length. Nuclear data is what licenses biparental ancestry, hybridization, recombination, inbreeding, and introgression claims.

---

## 8. Status and next step

**FINAL — APPROVED (v2.5).** All issues from six review rounds are resolved in-document: scope is fixed (4 checks + 6 figures, presentation changes folded in, geography/NeighborNet conditional); the metadata blocker is DUSE-only with sample-level provenance (§2.5); index-hopping and duplicate-library claims are corrected and framed as screens (§3.2); audit tolerances carry magnitude as well as rank criteria, separate baseline from sensitivity, and separate π/FST eligibility from haplotype-assignment eligibility (§6); likelihood mapping has a committed input, predeclared seeds/model/thresholds, and a within-organelle bootstrap-frequency NeighborNet trigger that is logically valid (§5); the RF and 1,000-replicate seeded downsampling procedures are fully specified (§6); coordinates are versioned under the supplementary run, never the frozen canonical file (§3.5); and the two documentation safeguards are Phase-1 gating deliverables (§6 items 5–6).

**Recommended immediate next step:** turn **Phase 1** (six gating deliverables — DUSE first, the identity audit, provenance re-verification, filtering sensitivity, the claim→analysis matrix, and the inheritance evidence) into an executable spec under `./specs/` with exact file paths, the §6 tolerances, the §2.5 sample-level provenance policy, and the §3.2 four-outcome identity rule — so execution does not relitigate this document.
