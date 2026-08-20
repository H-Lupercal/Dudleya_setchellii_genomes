# Downstream Sample QC Decisions

The Step 5 alignments remain available for audit, but the following samples are
ignored for downstream primary population-genetic analyses:

| Sample | cpDNA use | mtDNA use | Reason |
|---|---|---|---|
| `ABAB_MAD_LP_222_Du-589` | exclude | exclude | Extremely small FASTQ input and failed cpDNA/mtDNA coverage QC. |
| `CY_HUN_LP_265_Du-684` | exclude | exclude | Third-smallest FASTQ input in the run and failed mtDNA coverage QC; cpDNA coverage is usable, but this sample is excluded conservatively with the flagged low-input set. |
| `CY_RED_LP_202_Du-561` | exclude | exclude | Extremely small FASTQ input and failed cpDNA/mtDNA coverage QC. |

Downstream scripts for haploid variant calling, consensus FASTA generation,
all-sample alignments, PCA, phylogenetic trees, Fst, and
structure/admixture-style clustering should exclude these three samples from the
primary analysis set.

The machine-readable version of this decision is:

```text
dudleya_organelle_alignment_pipeline/results/06_all_sample_alignment/downstream_sample_qc_decisions.tsv
```
