# NC_085682.1 Chloroplast Comparison

This directory contains a focused comparison between the candidate *Dudleya setchellii*
chloroplast assembly and the public *Dudleya farinosa* chloroplast reference
`NC_085682.1`.

NCBI source record: https://www.ncbi.nlm.nih.gov/nuccore/NC_085682.1

## Purpose

The comparison answers four questions:

- Whole-plastome similarity and percent sequence divergence.
- Synonymous versus nonsynonymous single-nucleotide CDS substitutions.
- Gene content relative to `NC_085682.1`.
- Gene order and orientation relative to `NC_085682.1`.

## Inputs

- Candidate raw chloroplast FASTA:
  `../Dudleya_hifiasm_purged_manual_chloroplast.fa`
- Downloaded NCBI reference records:
  `NC_085682.1.fetched.fa` and `NC_085682.1.fetched.gb`
- Candidate draft annotation used only for annotation-count comparison:
  `../gateVER/chloroplast.best_nonredundant.annotation.tsv`

The script downloads current `NC_085682.1` FASTA and GenBank records from NCBI
EFetch. If network access is unavailable, the normal run fails instead of using
cached reference data. This is intentional so the committed comparison artifacts
come from NCBI, not from prior gateVER cache state.

## Dependencies

Install these on `PATH`:

- `python3`
- `blastn` from BLAST+
- `mafft`

No Python package dependencies are required.

## Reproduce The Analysis

Run from the repository root:

```bash
python3 nc085682_comparison/run_comparison.py
```

The script rewrites the comparison artifacts in this directory. It records the
NCBI Nuccore URL and EFetch URLs in `reference_metadata.tsv`.

For exploratory debugging only, the script has an `--allow-cached-fallback`
option. Do not use that option for the results committed in this repository.

## Method Summary

1. Read the raw candidate chloroplast sequence.
2. Trim the known terminal duplicate by keeping candidate positions `1..150274`.
   The removed duplicate suffix starts at position `150275`, matching the
   existing gateVER repeat evidence.
3. BLAST the raw candidate against `NC_085682.1` to reproduce the existing
   gateVER similarity result.
4. BLAST `NC_085682.1` against the de-duplicated candidate to identify the
   circular-origin rotation. In this run, reference position 1 maps to candidate
   position `140018`.
5. Rotate the de-duplicated candidate to the `NC_085682.1` origin.
6. Build a normalized BLAST projection using aligned `qseq` and `sseq` strings.
   This projection is the primary source for divergence, CDS substitution,
   gene-content, and gene-order calls.
7. Run MAFFT on the normalized pair as a diagnostic only. Large repeats and SSC
   orientation can make a linear whole-genome alignment overstate divergence, so
   MAFFT numbers are not used as final similarity calls.
8. Parse `NC_085682.1` GenBank features, project features onto the candidate,
   and classify CDS codon differences. Single-nucleotide codon changes are
   counted as synonymous or nonsynonymous. Multi-nucleotide codon changes are
   kept in `complex_codon_changes`.

## Key Outputs

- `similarity_results.tsv`: one-row summary of the similarity run.
- `summary.md`: narrative interpretation and headline numbers.
- `whole_genome_divergence.tsv`: raw BLAST and normalized projection divergence.
- `cds_substitutions.tsv`: CDS-level synonymous/nonsynonymous substitution table.
- `gene_content.tsv`: reference feature counts, projected candidate presence,
  and candidate draft-annotation counts.
- `gene_order.tsv`: projected gene order, orientation, and IR/SSC notes.
- `reference_metadata.tsv`: accession, length, source, fetch date, and NCBI URLs.
- `verification_metrics.tsv`: compact metrics used by the verification checks.
- `NC_085682.1.fetched.fa` and `NC_085682.1.fetched.gb`: the FASTA and
  GenBank records downloaded from NCBI for this comparison run.

## Current Result Snapshot

- Raw BLAST weighted identity: `99.359%`, about `0.641%` divergence.
- Normalized similarity excluding gaps: `99.556017%`, about `0.443983%`
  divergence.
- Normalized similarity counting reference deletions: `99.365529%`, about
  `0.634471%` divergence.
- Normalized similarity counting unmapped reference bases as differences:
  `99.193527%`, about `0.806473%` divergence.
- CDS clean substitution totals: `100` synonymous and `93` nonsynonymous
  single-nucleotide substitutions across `79` clean projected CDS rows.
- Gene content: all `226/226` `NC_085682.1` feature rows are present by
  projection.
- Gene order: `47` rows are shifted/reordered after circular normalization and
  `31` rows project on the opposite strand, consistent with SSC orientation and
  IR-copy effects rather than broad gene-content loss.

## Interpretation Notes

The raw candidate is `176,964 bp`, while `NC_085682.1` is `150,780 bp`. The
candidate includes an approximately `26,690 bp` terminal duplicate, so the
biologically meaningful comparison uses the de-duplicated `150,274 bp` candidate
representation.

The synonymous/nonsynonymous counts are codon-consequence counts from pairwise
projected CDS sequences. They are not model-based dN/dS or Ka/Ks estimates.

The current candidate annotation is a homology-transfer draft. `gene_content.tsv`
separates sequence presence by projection from differences in draft annotation
copy count.
