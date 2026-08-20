# Concatenated Organelle Bootstrap Design

## Objective

Build one combined organelle alignment by appending each sample's mtDNA
callable-consensus sequence to the end of the same sample's cpDNA
callable-consensus sequence. Prepare, but do not start, a 10,000-replicate
IQ-TREE analysis of that combined alignment. Preserve every existing script and
all existing cpDNA- and mtDNA-only outputs.

## Confirmed Inputs

- cpDNA alignment: 275 samples, 124,538 sites.
- mtDNA alignment: 275 samples, 44,930 sites.
- The sample identifiers and their order match exactly between alignments.
- Combined alignment: 275 samples, 169,468 sites per sample.
- Coordinate boundary: cpDNA 1-124,538; mtDNA 124,539-169,468.

## Analysis Semantics

Each biological sample remains one phylogenetic taxon. Its output record is:

```text
combined_sequence[sample] = cpDNA_sequence[sample] + mtDNA_sequence[sample]
```

The combined alignment will be analyzed as one unpartitioned DNA alignment
under `GTR+F+G4`. No partition file will be passed to IQ-TREE. IQ-TREE will use
10,000 ultrafast bootstrap replicates with BNNI correction and 14 CPU threads.
The cpDNA-first order is a coordinate convention; reversing the order would not
add an independent biological analysis.

## Code and File Isolation

Existing scripts and modules will not be edited. New files will provide the new
behavior:

1. A concatenation module that validates and combines the two FASTA alignments.
2. A thin concatenation command-line script.
3. A copied and renamed 10,000-bootstrap runner configured for the combined
   alignment and its own output directory.
4. Focused tests for concatenation and fixed runner settings.

The generic existing phylogenetic-tree implementation will be reused without
modification because its summary-driven input supports a single combined row.

## Concatenation Outputs

The concatenation stage will write to:

```text
dudleya_organelle_alignment_pipeline/results/22_concatenated_consensus/
```

It will contain:

- `cpDNA_mtDNA.primary.concatenated_consensus.fa`: the 275 combined records.
- `primary.callable_consensus_summary.tsv`: compatibility metadata consumed by
  the existing tree implementation.
- `primary.concatenated_consensus_summary.tsv`: explicit concatenation metrics.
- `primary.concatenated_consensus_report.md`: human-readable provenance and
  boundary information.

## Bootstrap Outputs

The prepared runner will write only to:

```text
dudleya_organelle_alignment_pipeline/results/23_concatenated_bootstrap_phylogenetic_tree_10000/
```

It will be configured with `--bootstrap-replicates 10000`, `--threads 14`, and
the Stage 22 consensus directory. This runner will not be executed until the
user explicitly gives the go-ahead after reviewing the completed
concatenation.

## Validation and Failure Handling

Concatenation must stop without writing a successful result if any of these
conditions holds:

- either alignment is empty;
- a sample identifier is duplicated;
- the sample identifier sets differ;
- records have inconsistent lengths within either alignment;
- any output sequence length differs from 169,468 sites;
- an output record is not the exact cpDNA sequence followed by the exact mtDNA
  sequence for that identifier.

The cpDNA record order will determine output order. Missing and ambiguous bases
will be retained unchanged. The script will report the sample count, component
lengths, combined length, coordinate boundary, and total missing-base count.

## Testing and Execution Boundary

Tests will be written before implementation and will cover:

- exact cpDNA-then-mtDNA concatenation by sample identifier;
- preservation of cpDNA order when mtDNA input order differs;
- rejection of mismatched and duplicate sample identifiers;
- rejection of inconsistent sequence lengths;
- generation of tree-runner arguments containing 10,000 replicates, 14
  threads, and isolated input/output directories.

After the tests pass, only the concatenation stage will run. Completion requires
fresh verification that the FASTA contains exactly 275 unique records, every
record is 169,468 sites, the seam is exact for every sample, and the reports
record the 1-124,538 and 124,539-169,468 boundaries. The bootstrap process must
not be present after this step.

## Interpretation Caveat

Because the analysis is deliberately unpartitioned, one substitution model and
rate distribution will describe both organelles. cpDNA supplies approximately
73.5% of the concatenated columns, although each organelle's effective
phylogenetic influence depends on its variable and informative sites. The final
tree represents combined organellar signal and may conceal cpDNA/mtDNA conflict;
the separate trees remain necessary companion results.
