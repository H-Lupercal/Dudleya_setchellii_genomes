# gateVER annotation integrity report

## Scope

This report uses only files inside `gateVER/`: input FASTAs, best nonredundant annotation TSV/GFF3 files, cached GenBank references, and cached whole-genome BLAST summaries. The older `dudleya_organelle_qc/` directory is intentionally ignored.

No existing gateVER FASTA, TSV, or GFF3 annotation file was modified; this directory contains report-only derived QC artifacts for validation, not publication-ready manual curation.

## Input and annotation summary

| Organelle | Contig | Length bp | GC % | gateVER CDS | tRNA | rRNA |
|---|---|---:|---:|---:|---:|---:|
| chloroplast | `ptg000216l_1` | 176,964 | 37.77 | 85 | 30 | 4 |
| mitochondria | `ptg000317l_1` | 243,359 | 45.517 | 34 | 13 | 2 |

## Validation procedures completed

| Category | Artifact | Status | Result |
|---|---|---|---|
| CDS integrity | `cds_integrity.tsv` | DONE | 119 CDS checked across cpDNA and mtDNA |
| Expected genes | `expected_gene_presence.tsv` | DONE | {'OK': 177, 'MISSING_EXPECTED': 4, 'CANDIDATE_ONLY_OR_RARE_REFERENCE': 5} |
| Mitochondrial core genes | `mitochondrial_core_gene_presence.tsv` | DONE | {'PRESENT': 23, 'NOT_ANNOTATED_IN_GATEVER': 11} |
| Chloroplast IR features | `chloroplast_ir_feature_check.tsv` | DONE | {'ONE_IR_COPY_ANNOTATED': 35, 'IR_BOUNDARY_FEATURE': 6} |
| tRNA/rRNA plausibility | `trna_rrna_plausibility.tsv` | DONE | {'OK': 57, 'CANDIDATE_ONLY_OR_RARE_REFERENCE': 1, 'MISSING_EXPECTED': 1} |
| Whole-genome synteny | `whole_genome_synteny_check.tsv` | DONE | {'HIGH_COLLINEARITY_OR_FEW_BLOCKS': 6, 'PARTIAL_OR_REARRANGED_REFERENCE_MATCH': 3, 'LOWER_SUPPORT_REFERENCE': 7} |
| Cross-organelle identity | `cross_organelle_validation.tsv` | DONE | {'PASS_ORGANELLE_LABEL_SUPPORTED': 2} |
| Repeat structure | `repeat_structure_check.tsv` | DONE | Top large repeats reported for chloroplast and mitochondria |
| Read-backed validation availability | `read_graph_tool_availability.tsv` | DATA_NOT_AVAILABLE | No read/graph files found in repository |
| Manual visual review targets | `manual_review_targets.tsv` | DONE | 37 targets listed |

## Main findings

- Chloroplast identity remains strong in gateVER: top cached Dudleya plastome hit covers 100.00% of the candidate at 99.565% weighted identity.
- Mitochondrial identity remains strong at family/organelle level in gateVER: top cached Crassulaceae mitogenome hit covers 70.99% of the candidate at 97.69% weighted identity.
- gateVER chloroplast input has a terminal direct duplicate of about 26,690 bp; the plastome structure table includes a deduplicated comparison row without modifying the input FASTA.
- No read alignment, raw read, or assembly graph files were found in the repository, so coverage, SNP/indel hotspot, and read-backed mtDNA repeat/junction checks cannot be performed from current repo contents.
- gateVER annotations are useful draft homology-transfer calls. The integrity tables validate plausibility and identify review targets; they do not manually curate gene models.

## CDS integrity summary

| Organelle | CDS checked | PASS | WARN | REVIEW |
|---|---:|---:|---:|---:|
| chloroplast | 85 | 51 | 25 | 9 |
| mitochondria | 34 | 11 | 14 | 9 |

`PASS` means the current continuous gateVER CDS interval has ATG, length divisible by 3, no internal stop, and a terminal stop. `WARN` means the call is homology-supported but has boundary, split-gene, protein-transfer, or plant mitochondrial caveats. `REVIEW` means the model should not be treated as a curated CDS without manual inspection.

## Additional validation summaries

- Expected gene presence status counts: {'OK': 177, 'MISSING_EXPECTED': 4, 'CANDIDATE_ONLY_OR_RARE_REFERENCE': 5}.
- Mitochondrial core gene status counts: {'PRESENT': 23, 'NOT_ANNOTATED_IN_GATEVER': 11}.
- Chloroplast IR feature status counts: {'ONE_IR_COPY_ANNOTATED': 35, 'IR_BOUNDARY_FEATURE': 6}.
- tRNA/rRNA plausibility status counts: {'OK': 57, 'CANDIDATE_ONLY_OR_RARE_REFERENCE': 1, 'MISSING_EXPECTED': 1}.
- Whole-genome synteny status counts: {'HIGH_COLLINEARITY_OR_FEW_BLOCKS': 6, 'PARTIAL_OR_REARRANGED_REFERENCE_MATCH': 3, 'LOWER_SUPPORT_REFERENCE': 7}.
- Cross-organelle validation status counts: {'PASS_ORGANELLE_LABEL_SUPPORTED': 2}.

## Repeat structure highlights

| Organelle | Rank | Length bp | Identity % | Coordinates | Orientation | Interpretation |
|---|---:|---:|---:|---|---|---|
| chloroplast | 1 | 26702 | 99.944 | 150275-176964 vs 1-26700 | forward | terminal redundancy candidate |
| chloroplast | 2 | 26702 | 99.944 | 1-26700 vs 150275-176964 | forward | terminal redundancy candidate |
| chloroplast | 3 | 25742 | 99.953 | 114282-140017 vs 71834-97569 | reverse | inverted repeat candidate |
| chloroplast | 4 | 25742 | 99.953 | 71834-97569 vs 114282-140017 | reverse | inverted repeat candidate |
| mitochondria | 1 | 69706 | 99.999 | 153662-223367 vs 84115-153820 | forward | large internal/direct repeat candidate |
| mitochondria | 2 | 69706 | 99.999 | 84115-153820 vs 153662-223367 | forward | large internal/direct repeat candidate |

## Gene content comparison

See `gene_content_comparison.tsv` and `expected_gene_presence.tsv` for full candidate-versus-reference counts. Reference feature-copy counts include IR-duplicated chloroplast genes where present; unique-gene ranges collapse duplicated names.

## Plastome structure

| Sample | Length bp | Terminal duplicate bp | IR bp | LSC bp | SSC bp |
|---|---:|---:|---:|---:|---:|
| D_setchellii_candidate_chloroplast_raw | 176,964 | 26690 | 25742 | 108780 | 16712 |
| D_setchellii_candidate_chloroplast_terminal_deduplicated | 150,274 | 0 | 25742 | 82090 | 16712 |
| PX244394.1 | 166,371 | 0 | 30993 | 90859 | 13526 |

## Generated files

- `cds_integrity.tsv`
- `annotation_flags.tsv`
- `expected_gene_presence.tsv`
- `mitochondrial_core_gene_presence.tsv`
- `chloroplast_ir_feature_check.tsv`
- `trna_rrna_plausibility.tsv`
- `whole_genome_synteny_check.tsv`
- `cross_organelle_validation.tsv`
- `repeat_structure_check.tsv`
- `read_graph_tool_availability.tsv`
- `manual_review_targets.tsv`
- `validation_checklist.tsv`
- `gene_content_comparison.tsv`
- `plastome_structure_comparison.tsv`
- `assembly_flags.tsv`
- `manifest.json`

## Recommended use

Use this package to support validation claims: organelle identity, broad assembly plausibility, expected gene content, cpDNA quadripartite structure, and annotation plausibility. For publication/submission-grade genomes, separate work would still be needed for read-backed mtDNA structure validation and manual correction of `WARN`/`REVIEW` CDS models.
