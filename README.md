# Dudleya setchellii Organelle Reference Verification

This repository contains verified *Dudleya setchellii* chloroplast and
mitochondrial reference assets, supporting evidence, downloaded sequencing data,
and a local copy of the Dudleya conservation-genomics pipeline.

The canonical organelle verification package is:

```text
dudleya_organelle_reference_verification/
```

Older analysis folders were merged into that package and removed to avoid
multiple competing sources of truth.

## Key Contents

| Path | Purpose |
|---|---|
| `dudleya_organelle_reference_verification/` | Canonical cpDNA/mtDNA references, annotations, identity evidence, NC_085682 chloroplast comparison, independent BLAST QC, and annotation-integrity checks. |
| `dudleya_organelle_alignment_pipeline/` | Reproducible cpDNA/mtDNA FASTQ-processing and population-genomics pipeline, including QC, alignments, variant calling, PCA, trees, admixture, and Fst outputs. |
| `dudleya_conservation_genomics_pipeline/` | Local credited copy of the SCU Dudleya conservation genomics pipeline originally published at `https://github.com/evanhackstadt/dudleya`. |
| `genomicsDrive_data_dump/` | Downloaded sequencing FASTQ data. |
| `ORGANELLE_POPGEN_WORK_PLAN.md` | Work plan for cpDNA/mtDNA population analysis from the FASTQ data. |
| `Dudleya_hifiasm_purged_manual_chloroplast.fa` | Original raw chloroplast candidate assembly retained for provenance. |
| `Dudleya_hifiasm_purged_manual_mitochondria.fa` | Original raw mitochondrial candidate assembly retained for provenance. |

## Final Organelle Population-Genomics Results

The professor's requested deliverables are complete for cpDNA and mtDNA. Start
with the handoff report, then use the specific output files below as needed.
Generated analysis products remain local and are intentionally ignored by git;
the three professor-facing result indexes under `results/` stay visible.

| File | What it is |
|---|---|
| [`dudleya_organelle_alignment_pipeline/results/PROFESSOR_HANDOFF.md`](dudleya_organelle_alignment_pipeline/results/PROFESSOR_HANDOFF.md) | Short professor-facing completion note with the final goal, what was produced, and caveats. |
| [`dudleya_organelle_alignment_pipeline/results/organelle_population_report.md`](dudleya_organelle_alignment_pipeline/results/organelle_population_report.md) | Integrated methods/results report covering references, QC, alignments, PCA, ML trees, admixture, Fst, and caveats. |
| [`dudleya_organelle_alignment_pipeline/results/final_deliverables_manifest.tsv`](dudleya_organelle_alignment_pipeline/results/final_deliverables_manifest.tsv) | Machine-readable list of the final deliverables and notes. |

### Alignments

| File | What it is |
|---|---|
| [`dudleya_organelle_alignment_pipeline/results/11_callable_consensus/cpDNA.primary.callable_consensus.fa`](dudleya_organelle_alignment_pipeline/results/11_callable_consensus/cpDNA.primary.callable_consensus.fa) | cpDNA full callable-site consensus alignment: 275 samples x 124,538 callable sites. |
| [`dudleya_organelle_alignment_pipeline/results/11_callable_consensus/mtDNA.primary.callable_consensus.fa`](dudleya_organelle_alignment_pipeline/results/11_callable_consensus/mtDNA.primary.callable_consensus.fa) | mtDNA high-confidence unique-track callable consensus alignment: 275 samples x 44,930 callable sites. |
| [`dudleya_organelle_alignment_pipeline/results/10_snp_alignment/cpDNA.primary.snp_alignment.fa`](dudleya_organelle_alignment_pipeline/results/10_snp_alignment/cpDNA.primary.snp_alignment.fa) | cpDNA filtered haploid SNP alignment with 2,015 SNP sites. |
| [`dudleya_organelle_alignment_pipeline/results/10_snp_alignment/mtDNA.primary.snp_alignment.fa`](dudleya_organelle_alignment_pipeline/results/10_snp_alignment/mtDNA.primary.snp_alignment.fa) | mtDNA filtered haploid SNP alignment with 146 SNP sites. |

### PCA

| File | What it is |
|---|---|
| [`dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.png`](dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.png) | cpDNA PCA plot; PC1 explains 36.62 percent and PC2 explains 14.65 percent. |
| [`dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.png`](dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.png) | mtDNA PCA plot; PC1 explains 34.48 percent and PC2 explains 14.06 percent. |
| [`dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.coordinates.tsv`](dudleya_organelle_alignment_pipeline/results/15_pca/cpDNA.primary.pca.coordinates.tsv) | cpDNA sample coordinates and metadata for replotting or labeling PCA results. |
| [`dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.coordinates.tsv`](dudleya_organelle_alignment_pipeline/results/15_pca/mtDNA.primary.pca.coordinates.tsv) | mtDNA sample coordinates and metadata for replotting or labeling PCA results. |

### Phylogenetic Trees

| File | What it is |
|---|---|
| [`dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile`](dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile) | cpDNA IQ-TREE maximum-likelihood tree, GTR+F+G4, 1,000 ultrafast bootstraps with BNNI. |
| [`dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile`](dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile) | mtDNA IQ-TREE maximum-likelihood tree, GTR+F+G4, 1,000 ultrafast bootstraps with BNNI. |
| [`dudleya_organelle_alignment_pipeline/results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png`](dudleya_organelle_alignment_pipeline/results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png) | Rendered cpDNA ML tree figure with bootstrap support. |
| [`dudleya_organelle_alignment_pipeline/results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png`](dudleya_organelle_alignment_pipeline/results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png) | Rendered mtDNA ML tree figure with bootstrap support. |

### Admixture And Population Statistics

| File | What it is |
|---|---|
| [`dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png`](dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png) | cpDNA ADMIXTURE-style structure plot; replicate-based best K=8. |
| [`dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png`](dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png) | mtDNA ADMIXTURE-style structure plot; replicate-based best K=8. |
| [`dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/primary.admixture_report.md`](dudleya_organelle_alignment_pipeline/results/18_admixture_replicates/primary.admixture_report.md) | K=1..8, five-replicate CV-error sweep and haploid pseudo-diploid encoding note. |
| [`dudleya_organelle_alignment_pipeline/results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv`](dudleya_organelle_alignment_pipeline/results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv) | cpDNA pairwise Fst table: 561 comparisons across 34 resolved populations. |
| [`dudleya_organelle_alignment_pipeline/results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv`](dudleya_organelle_alignment_pipeline/results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv) | mtDNA pairwise Fst table: 561 comparisons across 34 resolved populations. |
| [`dudleya_organelle_alignment_pipeline/results/17_population_genetics/cpDNA.primary.population_genetics.population_summary.tsv`](dudleya_organelle_alignment_pipeline/results/17_population_genetics/cpDNA.primary.population_genetics.population_summary.tsv) | cpDNA population summary with sample counts, haplotypes, diversity, nucleotide diversity, and private variants. |
| [`dudleya_organelle_alignment_pipeline/results/17_population_genetics/mtDNA.primary.population_genetics.population_summary.tsv`](dudleya_organelle_alignment_pipeline/results/17_population_genetics/mtDNA.primary.population_genetics.population_summary.tsv) | mtDNA population summary with sample counts, haplotypes, diversity, nucleotide diversity, and private variants. |
| [`dudleya_organelle_alignment_pipeline/results/13_tool_audit/primary.tool_audit_report.md`](dudleya_organelle_alignment_pipeline/results/13_tool_audit/primary.tool_audit_report.md) | Installed/verified bioinformatics and visualization tool audit. |

Main interpretation caveats: trees are ML rather than a standalone rendered NJ
tree; Fst is a custom haploid Nei-style estimate over informative SNPs; and
mtDNA analyses use the 44,930 bp high-confidence unique track, leaving
repeat-rich mtDNA regions out of population-genetic interpretation.

## Canonical References

Use these files for the next mapping and population-analysis steps:

```text
dudleya_organelle_reference_verification/references/chloroplast.normalized.fa
dudleya_organelle_reference_verification/references/mitochondria.fa
dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa
```

The chloroplast reference has been terminal-deduplicated and rotated to the
`NC_085682.1` origin. The mitochondrial reference keeps the verified
mitochondrial candidate sequence with a stable `mitochondria` FASTA header.

## Canonical Annotations

```text
dudleya_organelle_reference_verification/annotations/chloroplast.gff3
dudleya_organelle_reference_verification/annotations/chloroplast.annotation.tsv
dudleya_organelle_reference_verification/annotations/mitochondria.gff3
dudleya_organelle_reference_verification/annotations/mitochondria.annotation.tsv
```

These are homology-transfer draft annotations. They are useful for feature-aware
interpretation and QC, but they are not curated GenBank-submission annotations.

## Evidence Summary

The chloroplast assembly is strongly supported as a Dudleya chloroplast genome.
Whole-genome comparisons against complete Dudleya plastomes show about
99.3-99.6 percent weighted nucleotide identity. The `NC_085682.1` comparison
also shows that the raw chloroplast candidate contains a terminal duplicate; the
normalized reference in `references/chloroplast.normalized.fa` removes that
duplicate and rotates the circular genome to the public reference origin.

The mitochondrial assembly is supported as a Crassulaceae mitochondrial genome.
Whole-genome comparisons against related mitochondrial references show about
97.3-97.7 percent weighted nucleotide identity across large portions of the
candidate. Cross-organelle checks found no evidence that the chloroplast and
mitochondrial labels are swapped.

The annotation-integrity checks classify the best nonredundant CDS calls as:

```text
chloroplast: 51 PASS, 25 WARN, 9 REVIEW
mitochondria: 11 PASS, 14 WARN, 9 REVIEW
```

See:

```text
dudleya_organelle_reference_verification/README.md
dudleya_organelle_reference_verification/annotation_integrity_checks/report.md
```

## Pipeline Status

The cpDNA/mtDNA FASTQ-processing and downstream population-genomics pipeline is
implemented here:

```text
dudleya_organelle_alignment_pipeline/
```

The final professor-facing deliverables have already been run through alignment,
haploid variant calling, SNP/callable alignments, PCA, maximum-likelihood trees,
admixture-style structure plots, and population-genetic summaries. Recreate or
validate the earlier manifest/reference/pilot/QC and smoke variant-calling steps
from the repo root with:

```bash
python3 -m unittest \
  dudleya_organelle_alignment_pipeline.tests.test_manifest \
  dudleya_organelle_alignment_pipeline.tests.test_prepare_reference_and_pilot \
  dudleya_organelle_alignment_pipeline.tests.test_pilot_alignment \
  dudleya_organelle_alignment_pipeline.tests.test_analysis_masks \
  dudleya_organelle_alignment_pipeline.tests.test_all_sample_alignment \
  dudleya_organelle_alignment_pipeline.tests.test_downstream_sample_set \
  dudleya_organelle_alignment_pipeline.tests.test_variant_calling -v
python3 dudleya_organelle_alignment_pipeline/scripts/build_sample_manifest.py
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/prepare_reference_and_pilot.py
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_pilot_alignment.py
python3 dudleya_organelle_alignment_pipeline/scripts/build_analysis_masks.py
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_all_sample_alignment.py
python3 dudleya_organelle_alignment_pipeline/scripts/build_downstream_sample_set.py
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_variant_calling.py \
  --run-label smoke \
  --sample-id ABAB_MAD_LP_223_Du-590 \
  --sample-id ABAB_MAD_LP_322_Du-593 \
  --sample-id ABAB_MAD_LP_323_Du-594 \
  --sample-id ABAB_MAD_LP_324_Du-595 \
  --sample-id ABAB_MAD_LP_325_Du-596
```

Use this table for the primary paired-end cpDNA/mtDNA alignment:

```text
dudleya_organelle_alignment_pipeline/results/00_manifest/analysis_samples.tsv
```

Two manually verified missing-mate samples are excluded from that primary table
and documented in:

```text
dudleya_organelle_alignment_pipeline/results/00_manifest/excluded_samples.tsv
```

The current pilot table is:

```text
dudleya_organelle_alignment_pipeline/results/01_reference_pilot/pilot_samples.tsv
```

A local repo-specific bioinformatics environment has been created at
`.tools/bioconda-env/`. It is intentionally ignored by git; recreate it from
`dudleya_organelle_alignment_pipeline/environment.yml` if needed. Step 2 now
finds `bwa`, `samtools`, `fastp`, `fastqc`, `multiqc`, and `bcftools` when run
with that environment on `PATH`, and it has created the `samtools faidx` and
`bwa index` files for `dudleya_cp_mt.fa`.

Step 3 pilot alignment writes filtered organelle BAMs and depth files under:

```text
dudleya_organelle_alignment_pipeline/results/02_pilot_alignment/
```

Only the small top-level pilot summaries are intended to be kept in git.
The current pilot run summarized 15 samples and 30 sample-by-organelle rows.
After correcting the `samtools depth` quality flags, median breadth at `>=1x`
is about `1.000` for chloroplast and `0.960` for mitochondria. The remaining
mtDNA concern is not broad absence of coverage; it is that repeat-rich mtDNA
regions often have low mapping quality and need a separate unique-placement
mask before variant calling.

The focused mtDNA investigation report is:

```text
dudleya_organelle_alignment_pipeline/results/03_mtdna_investigation/mtdna_investigation_report.md
```

The focused cpDNA verification report is:

```text
dudleya_organelle_alignment_pipeline/results/04_cpdna_investigation/cpdna_verification_report.md
```

The cpDNA pilot verification supports moving forward with all-sample chloroplast
processing. The main cpDNA caution is standard chloroplast inverted-repeat
handling: the normalized reference has a 25,742 bp reverse repeat pair at
`82091-107826` and `124539-150274`, so downstream SNP analyses should mask one
IR copy or otherwise avoid counting duplicated IR sequence twice.

Step 4 defines those analysis rules as machine-readable BED/TSV files:

```text
dudleya_organelle_alignment_pipeline/results/05_analysis_masks/
```

Key Step 4 decisions:

- cpDNA sample QC can use the full chloroplast reference.
- cpDNA population-genetic outputs should use
  `cpdna_population_sites.bed`, which keeps one IR copy and excludes the
  duplicate IR copy `124539-150274`.
- mtDNA sample QC should use the whole-reference permissive coverage track.
- mtDNA variant calling and population-genetic outputs should use
  `mtdna_high_confidence_unique_regions.bed`, currently 44,930 bp supported by
  high-MAPQ pilot evidence.

Step 5 writes the all-sample alignment and track-aware QC summaries under:

```text
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/
```

Step 6 writes the primary downstream include/exclude sample set under:

```text
dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/
```

Step 7 writes raw haploid cpDNA/mtDNA variant-calling outputs under:

```text
dudleya_organelle_alignment_pipeline/results/08_variant_calling/
```

## Limitations

The verification evidence is FASTA-based. Read-backed coverage, SNP/indel
review, and mitochondrial repeat/junction validation require mapping the
downloaded FASTQ data back to the canonical references.

No license has been declared for these data. Contact the project owners before
reuse or redistribution beyond normal review and citation.
