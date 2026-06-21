# NC_085682.1 Chloroplast Comparison

## Reference and Candidate
- Reference: NC_085682.1, Dudleya farinosa chloroplast, complete genome. (150780 bp; LOCUS date 01-FEB-2024).
- Reference source: live NCBI EFetch; fetch date recorded as 2026-06-21.
- Candidate raw plastome: 176964 bp; ambiguous bases 0.
- Candidate primary comparison plastome: terminal de-duplicated to 150274 bp by trimming the suffix starting at position 150275.
- Rotation/orientation: BLAST anchor at candidate position 140018 in + orientation, pident 97.978, length 7664 bp.

## Whole-Genome Divergence
- Raw BLAST compatibility check: query coverage 0.999486, reference coverage 0.998269, weighted identity 99.359%, HSPs 9.
- One-row similarity result file: `similarity_results.tsv`.
- Normalized BLAST projection: 150519 mapped reference bp, 261 unmapped reference bp, 667 mismatches, and 288 candidate deletion bases versus the reference.
- Identity excluding gaps: 99.556017%; divergence excluding gaps: 0.443983%.
- Identity counting reference deletions: 99.365529%; divergence counting reference deletions: 0.634471%.
- Identity counting unmapped reference bases as differences: 99.193527%; divergence on that stricter denominator: 0.806473%.
- Diagnostic MAFFT linear alignment, not used for final calls because repeats/SSC orientation can mislead it: 153844 columns, 6976 mismatches, 6634 indel columns.

## CDS Substitutions
- Shared/projected CDS rows: 84.
- CDS included in clean single-codon substitution totals: 79.
- Clean single-nt synonymous substitutions: 100.
- Clean single-nt nonsynonymous substitutions: 93.
- Complex codon changes kept separate: 3; amino-acid-changing codons across all projected CDS rows: 119.
- Reference CDS translation validation: 84 matches and 0 mismatches among 84 CDS with /translation qualifiers.

## Gene Content and Order
- Gene-content rows with full projection support: 226/226.
- Rows present by alignment but differing from the current draft annotation count: 37. This is expected in IR and draft-transfer edge cases.
- Gene-order rows flagged as shifted/reordered after circular normalization: 47.
- Gene-order rows projected on the opposite strand: 31, consistent with an SSC-orientation difference and IR-copy effects. Inspect `gene_order.tsv` for exact rows.

## Interpretation
The candidate chloroplast is very close to NC_085682.1 after removing the terminal duplicate and rotating to the same circular origin. The raw BLAST result remains consistent with existing gateVER evidence, while the normalized BLAST projection gives the biologically more meaningful divergence estimate. Gene content is mostly present by alignment; differences from the draft annotation should be interpreted as annotation completeness/copy-number issues unless `gene_content.tsv` marks projection absence. Synonymous/nonsynonymous totals are codon-consequence counts from projected shared CDS and are not model-based dN/dS or Ka/Ks estimates.
