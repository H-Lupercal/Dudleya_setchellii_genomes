# Dudleya Organelle Reference Verification

This package is the canonical chloroplast and mitochondrial reference package
for the Dudleya organelle/population-analysis work. It consolidates the former
tool-output and comparison folders into one descriptive layout.

## Contents

| Path | Purpose |
|---|---|
| `references/` | Mapping-ready chloroplast, mitochondrial, and combined cp/mt FASTA references. |
| `annotations/` | Best nonredundant draft GFF3 and TSV annotations for cpDNA and mtDNA. |
| `evidence/identity/` | Whole-genome identity and cross-organelle summary tables. |
| `evidence/independent_blast_qc/` | Independent BLAST QC report, whole-genome tables, and marker-hit tables. |
| `evidence/nc085682_chloroplast_comparison/` | Focused comparison between the Dudleya chloroplast candidate and public `NC_085682.1`. |
| `annotation_integrity_checks/` | Annotation and assembly plausibility checks, review targets, and validation tables. |

## Mapping References

Use these for FASTQ mapping:

```text
references/chloroplast.normalized.fa
references/mitochondria.fa
references/dudleya_cp_mt.fa
```

`chloroplast.normalized.fa` is terminal-deduplicated and rotated to the
`NC_085682.1` origin. `mitochondria.fa` is the verified mitochondrial candidate
with a stable `mitochondria` FASTA header. `dudleya_cp_mt.fa` combines both
records for first-pass organelle read screening.

## Draft Annotations

```text
annotations/chloroplast.gff3
annotations/chloroplast.annotation.tsv
annotations/mitochondria.gff3
annotations/mitochondria.annotation.tsv
```

These annotations are homology-transfer drafts. Treat them as feature evidence
and review guides, not as final curated GenBank-submission annotations.

## Evidence Summary

The chloroplast FASTA is strongly supported as a Dudleya chloroplast genome. It
covers nearly all complete Dudleya chloroplast references at about 99.3-99.6
percent weighted nucleotide identity.

The mitochondrial FASTA is supported as a Crassulaceae mitochondrial genome. It
covers large portions of related mitochondrial references at about 97.3-97.7
percent weighted nucleotide identity.

Cross-organelle comparisons do not support the chloroplast and mitochondrial
labels being swapped.

The focused `NC_085682.1` comparison shows that the raw chloroplast candidate has
an approximately 26.7 kb terminal duplicate. The normalized chloroplast FASTA in
`references/chloroplast.normalized.fa` removes that duplicate and rotates the
sequence to the public reference origin.

## Annotation Integrity

Current best nonredundant CDS classifications:

```text
chloroplast: 51 PASS, 25 WARN, 9 REVIEW
mitochondria: 11 PASS, 14 WARN, 9 REVIEW
```

See:

```text
annotation_integrity_checks/report.md
annotation_integrity_checks/validation_checklist.tsv
annotation_integrity_checks/manual_review_targets.tsv
```

`PASS` means the current continuous interval has a clean ATG-to-stop ORF with no
internal stop. `WARN` and `REVIEW` calls remain useful homology evidence, but
need boundary, split-gene, RNA-editing, or overlap review before curated
annotation use.

## Limitations

This package is based on FASTA and homology evidence. It does not yet contain
read-backed coverage, SNP/indel, or mitochondrial repeat-junction validation.
Those checks require mapping the downloaded FASTQ R1/R2 files to
`references/dudleya_cp_mt.fa`.
