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
| `dudleya_organelle_alignment_pipeline/` | Reproducible cpDNA/mtDNA FASTQ-processing and population-genomics pipeline, including QC, alignments, variant calling, PCA, trees, admixture, and Fst outputs. See `PROCESS.md` for the ordered stage index. |
| `initial_pipeline_run/` | Organized entry point for the initial pipeline output snapshot, with real copied final results and a provenance symlink back to `dudleya_organelle_alignment_pipeline/results/`. |
| `full_pipeline_run/` | Organized second full pipeline rerun using all 16 CPU threads, with run report, logs, metadata, and final result files. |
| `dudleya_conservation_genomics_pipeline/` | Local credited copy of the SCU Dudleya conservation genomics pipeline originally published at `https://github.com/evanhackstadt/dudleya`. |
| `genomicsDrive_data_dump/` | Downloaded sequencing FASTQ data. |
| `ORGANELLE_POPGEN_WORK_PLAN.md` | Work plan for cpDNA/mtDNA population analysis from the FASTQ data. |
| `Dudleya_hifiasm_purged_manual_chloroplast.fa` | Original raw chloroplast candidate assembly retained for provenance. |
| `Dudleya_hifiasm_purged_manual_mitochondria.fa` | Original raw mitochondrial candidate assembly retained for provenance. |

## Final Organelle Population-Genomics Results

There are two organized output folders:

| Run | Folder | Notes |
|---|---|---|
| Initial run | [`initial_pipeline_run/`](initial_pipeline_run/) | Entry point for the first result snapshot, with copied final outputs under `results/` and a provenance symlink back to the original pipeline results. |
| Second full rerun | [`full_pipeline_run/`](full_pipeline_run/) | Complete rerun using all 16 CPU threads; stages `00` through `20` completed and produced a 34G organized result folder. |

Use the second full rerun for current interpretation unless you are explicitly
comparing against the initial run.

### Challenge Checklist

This repository explicitly addresses the requested organelle population-genomics
challenge as follows:

| Challenge requirement | Primary result |
|---|---|
| Create an alignment of all samples for cpDNA | [`full_pipeline_run/results/11_callable_consensus/cpDNA.primary.callable_consensus.fa`](full_pipeline_run/results/11_callable_consensus/cpDNA.primary.callable_consensus.fa) and [`full_pipeline_run/results/10_snp_alignment/cpDNA.primary.snp_alignment.fa`](full_pipeline_run/results/10_snp_alignment/cpDNA.primary.snp_alignment.fa). |
| Create an alignment of all samples for mtDNA | [`full_pipeline_run/results/11_callable_consensus/mtDNA.primary.callable_consensus.fa`](full_pipeline_run/results/11_callable_consensus/mtDNA.primary.callable_consensus.fa) and [`full_pipeline_run/results/10_snp_alignment/mtDNA.primary.snp_alignment.fa`](full_pipeline_run/results/10_snp_alignment/mtDNA.primary.snp_alignment.fa). |
| Use the annotated cpDNA and mtDNA references | Mapping and analysis use the canonical references in [`dudleya_organelle_reference_verification/references/`](dudleya_organelle_reference_verification/references/) with annotations in [`dudleya_organelle_reference_verification/annotations/`](dudleya_organelle_reference_verification/annotations/). |
| PCA visualization | [`full_pipeline_run/results/15_pca/cpDNA.primary.pca.png`](full_pipeline_run/results/15_pca/cpDNA.primary.pca.png), [`full_pipeline_run/results/15_pca/mtDNA.primary.pca.png`](full_pipeline_run/results/15_pca/mtDNA.primary.pca.png), and their coordinate tables. |
| Phylogenetic tree | Maximum-likelihood IQ-TREE outputs in [`full_pipeline_run/results/19_bootstrap_phylogenetic_tree/`](full_pipeline_run/results/19_bootstrap_phylogenetic_tree/) and rendered figures in [`full_pipeline_run/results/20_bootstrap_tree_visualization/`](full_pipeline_run/results/20_bootstrap_tree_visualization/). |
| ADMIXTURE / structure plot with empirical K | Five-replicate K=1..8 runs in [`full_pipeline_run/results/18_admixture_replicates/`](full_pipeline_run/results/18_admixture_replicates/); best K=8 for both cpDNA and mtDNA by mean CV error. |
| Fst and population-genetic parameters | Pairwise Fst and population summaries in [`full_pipeline_run/results/17_population_genetics/`](full_pipeline_run/results/17_population_genetics/). |
| QC and Evan pipeline context | The pipeline includes manifest checks, pilot alignment, all-sample alignment QC, downstream sample filtering, tool audit, and stage logs; the local credited copy of Evan's Dudleya pipeline is retained in [`dudleya_conservation_genomics_pipeline/`](dudleya_conservation_genomics_pipeline/) for context. |

| File | What it is |
|---|---|
| [`full_pipeline_run/README.md`](full_pipeline_run/README.md) | Result index with links to the organized rerun outputs. |
| [`full_pipeline_run/FULL_PIPELINE_RUN_REPORT.md`](full_pipeline_run/FULL_PIPELINE_RUN_REPORT.md) | Narrative report covering what ran, issues encountered, fixes, outputs, and verification. |
| [`full_pipeline_run/logs/stage_status.tsv`](full_pipeline_run/logs/stage_status.tsv) | Stage ledger showing completion of `00_manifest` through `20_bootstrap_tree_visualization`. |
| [`full_pipeline_run/run_metadata.txt`](full_pipeline_run/run_metadata.txt) | Run start/finish times, thread count, Python environment, and result size. |

### Alignments

| File | What it is |
|---|---|
| [`full_pipeline_run/results/11_callable_consensus/cpDNA.primary.callable_consensus.fa`](full_pipeline_run/results/11_callable_consensus/cpDNA.primary.callable_consensus.fa) | cpDNA callable-site consensus alignment: 278 samples x 124,538 callable sites. |
| [`full_pipeline_run/results/11_callable_consensus/mtDNA.primary.callable_consensus.fa`](full_pipeline_run/results/11_callable_consensus/mtDNA.primary.callable_consensus.fa) | mtDNA high-confidence unique-track callable consensus alignment: 278 samples x 44,930 callable sites. |
| [`full_pipeline_run/results/10_snp_alignment/cpDNA.primary.snp_alignment.fa`](full_pipeline_run/results/10_snp_alignment/cpDNA.primary.snp_alignment.fa) | cpDNA filtered haploid SNP alignment with 2,022 SNP sites. |
| [`full_pipeline_run/results/10_snp_alignment/mtDNA.primary.snp_alignment.fa`](full_pipeline_run/results/10_snp_alignment/mtDNA.primary.snp_alignment.fa) | mtDNA filtered haploid SNP alignment with 146 SNP sites. |

### PCA

| File | What it is |
|---|---|
| [`full_pipeline_run/results/15_pca/cpDNA.primary.pca.png`](full_pipeline_run/results/15_pca/cpDNA.primary.pca.png) | cpDNA PCA plot; PC1 explains 37.04 percent and PC2 explains 14.45 percent. |
| [`full_pipeline_run/results/15_pca/mtDNA.primary.pca.png`](full_pipeline_run/results/15_pca/mtDNA.primary.pca.png) | mtDNA PCA plot; PC1 explains 34.43 percent and PC2 explains 14.03 percent. |
| [`full_pipeline_run/results/15_pca/cpDNA.primary.pca.coordinates.tsv`](full_pipeline_run/results/15_pca/cpDNA.primary.pca.coordinates.tsv) | cpDNA sample coordinates and metadata for replotting or labeling PCA results. |
| [`full_pipeline_run/results/15_pca/mtDNA.primary.pca.coordinates.tsv`](full_pipeline_run/results/15_pca/mtDNA.primary.pca.coordinates.tsv) | mtDNA sample coordinates and metadata for replotting or labeling PCA results. |

### Phylogenetic Trees

| File | What it is |
|---|---|
| [`full_pipeline_run/results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile`](full_pipeline_run/results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile) | cpDNA IQ-TREE maximum-likelihood tree, GTR+F+G4, 1,000 ultrafast bootstraps with BNNI. |
| [`full_pipeline_run/results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile`](full_pipeline_run/results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile) | mtDNA IQ-TREE maximum-likelihood tree, GTR+F+G4, 1,000 ultrafast bootstraps with BNNI. |
| [`full_pipeline_run/results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png`](full_pipeline_run/results/20_bootstrap_tree_visualization/cpDNA.primary.iqtree_ml_tree.png) | Rendered cpDNA ML tree figure with bootstrap support. |
| [`full_pipeline_run/results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png`](full_pipeline_run/results/20_bootstrap_tree_visualization/mtDNA.primary.iqtree_ml_tree.png) | Rendered mtDNA ML tree figure with bootstrap support. |

### ADMIXTURE And Population Statistics

| File | What it is |
|---|---|
| [`full_pipeline_run/results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png`](full_pipeline_run/results/18_admixture_replicates/cpDNA.primary.bestK8.structure.png) | cpDNA ADMIXTURE-style structure plot; replicate-based best K=8. |
| [`full_pipeline_run/results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png`](full_pipeline_run/results/18_admixture_replicates/mtDNA.primary.bestK8.structure.png) | mtDNA ADMIXTURE-style structure plot; replicate-based best K=8. |
| [`full_pipeline_run/results/18_admixture_replicates/primary.admixture_report.md`](full_pipeline_run/results/18_admixture_replicates/primary.admixture_report.md) | K=1..8, five-replicate CV-error sweep and haploid pseudo-diploid encoding note. |
| [`full_pipeline_run/results/18_admixture_replicates/mtDNA.primary.pseudo_diploid.excluded_samples.tsv`](full_pipeline_run/results/18_admixture_replicates/mtDNA.primary.pseudo_diploid.excluded_samples.tsv) | ADMIXTURE-only mtDNA exclusion record for the single all-missing SNP-genotype sample. |
| [`full_pipeline_run/results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv`](full_pipeline_run/results/17_population_genetics/cpDNA.primary.population_genetics.pairwise_fst.tsv) | cpDNA pairwise Fst table: 595 comparisons across 35 resolved populations. |
| [`full_pipeline_run/results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv`](full_pipeline_run/results/17_population_genetics/mtDNA.primary.population_genetics.pairwise_fst.tsv) | mtDNA pairwise Fst table: 595 comparisons across 35 resolved populations. |
| [`full_pipeline_run/results/17_population_genetics/cpDNA.primary.population_genetics.population_summary.tsv`](full_pipeline_run/results/17_population_genetics/cpDNA.primary.population_genetics.population_summary.tsv) | cpDNA population summary with sample counts, haplotypes, diversity, nucleotide diversity, and private variants. |
| [`full_pipeline_run/results/17_population_genetics/mtDNA.primary.population_genetics.population_summary.tsv`](full_pipeline_run/results/17_population_genetics/mtDNA.primary.population_genetics.population_summary.tsv) | mtDNA population summary with sample counts, haplotypes, diversity, nucleotide diversity, and private variants. |
| [`full_pipeline_run/results/13_tool_audit/primary.tool_audit_report.md`](full_pipeline_run/results/13_tool_audit/primary.tool_audit_report.md) | Installed/verified bioinformatics and visualization tool audit from the rerun environment. |

Main interpretation caveats: trees are ML rather than a standalone rendered NJ
tree; Fst is a custom haploid Nei-style estimate over informative SNPs;
ADMIXTURE is a haploid organelle clustering visualization using pseudo-diploid
encoding; and mtDNA analyses use the 44,930 bp high-confidence unique track,
leaving repeat-rich mtDNA regions out of population-genetic interpretation.

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

The cpDNA/mtDNA FASTQ-processing and population-genomics pipeline is implemented
in `dudleya_organelle_alignment_pipeline/` and has been run to completion across
all stages, from `results/00_manifest/` through
`results/20_bootstrap_tree_visualization/`.

- Ordered stage index (authoritative): [`dudleya_organelle_alignment_pipeline/PROCESS.md`](dudleya_organelle_alignment_pipeline/PROCESS.md).
- Per-stage usage and commands: [`dudleya_organelle_alignment_pipeline/README.md`](dudleya_organelle_alignment_pipeline/README.md).

To recreate or validate the upstream manifest/reference/pilot/QC and smoke
variant-calling stages from the repository root:

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

Authoritative sample tables:

```text
dudleya_organelle_alignment_pipeline/results/00_manifest/analysis_samples.tsv
dudleya_organelle_alignment_pipeline/results/00_manifest/excluded_samples.tsv
dudleya_organelle_alignment_pipeline/results/07_downstream_sample_set/included_samples.tsv
```

`analysis_samples.tsv` is the primary paired-end input; `excluded_samples.tsv`
records the two missing-mate exclusions; `included_samples.tsv` is the 275-sample
downstream set used for variant calling and all population-genetic analyses.

The local bioinformatics environment lives at `.tools/bioconda-env/` (git-ignored;
recreate from `dudleya_organelle_alignment_pipeline/environment.yml`). Two analysis
cautions carry through to the population-genetic tracks: the chloroplast inverted
repeat (a 25,742 bp pair at `82091-107826` and `124539-150274`) is masked to one
copy for SNP analyses, and mtDNA population genetics is restricted to the 44,930 bp
high-confidence unique track. The supporting pilot investigations are in
`results/03_mtdna_investigation/` and `results/04_cpdna_investigation/`.

## Limitations

The verification evidence is FASTA-based. Read-backed coverage, SNP/indel
review, and mitochondrial repeat/junction validation require mapping the
downloaded FASTQ data back to the canonical references.

No license has been declared for these data. Contact the project owners before
reuse or redistribution beyond normal review and citation.
