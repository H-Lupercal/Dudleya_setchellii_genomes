# Dudleya setchellii Organelle Genome Verification

This repository contains two candidate organelle FASTA assemblies for *Dudleya setchellii*, independent QC evidence, draft homology-transfer annotations, gateVER annotation-integrity checks, and IGV-ready visualization files. It is meant to make the verification process reproducible by hand. It is not a packaged workflow or a one-command pipeline.

## Repository Contents

The root FASTA files are `Dudleya_hifiasm_purged_manual_chloroplast.fa` and `Dudleya_hifiasm_purged_manual_mitochondria.fa`. The `dudleya_organelle_qc/` directory contains an independent BLAST-based QC and identity report. The `gateVER/` directory contains reference search summaries, whole-genome comparison summaries, proof JSON, draft annotation TSV and GFF3 files, a detailed report, and the `gateVER/annotation_integrity/` report package. The `nc085682_comparison/` directory contains a focused, reproducible comparison between the candidate *Dudleya setchellii* chloroplast assembly and public *Dudleya farinosa* chloroplast reference `NC_085682.1`, including sequence divergence, synonymous/nonsynonymous CDS substitutions, gene content, gene order, and rerun documentation. The `igv/` directory contains combined and per-organelle FASTA, GFF3, and IGV JSON files for visual inspection.

## Requirements

To reproduce the analysis manually, install BLAST+ with `makeblastdb`, `blastn`, and `tblastn`. You also need a way to fetch NCBI nucleotide and GenBank records, such as Entrez Direct, `curl` against NCBI E-utilities, or the NCBI web interface. IGV Desktop is optional, but useful for inspecting the draft annotations.

## Manual Replication

Start with the two root FASTA files and record each contig name, length, GC percentage, and ambiguous base count. Use the NCBI searches and accessions in `gateVER/reference_search_summary.json`, `gateVER/report.md`, and `dudleya_organelle_qc/report.md` to download complete Dudleya chloroplast references and complete Crassulaceae mitochondrial references. Build nucleotide BLAST databases for the query and reference FASTAs with `makeblastdb`.

Run four whole-genome `blastn` comparisons: the chloroplast query against Dudleya chloroplast references, the chloroplast query against Crassulaceae mitochondrial references, the mitochondrial query against Crassulaceae mitochondrial references, and the mitochondrial query against Dudleya chloroplast references. Summarize each comparison by union query coverage, union reference coverage, weighted percent identity, HSP count, and aligned bases. Compare those summaries with the committed TSV and JSON evidence in `dudleya_organelle_qc/` and `gateVER/`.

Then extract CDS nucleotide and protein markers from the closest available reference records, especially `PX244389.1` for chloroplast markers and `PV256627.1` for mitochondrial markers. Align nucleotide markers with `blastn` and protein markers with `tblastn`, then compare marker support with `dudleya_organelle_qc/report.md` and the annotation evidence in `gateVER/`.

For draft annotation, transfer homologous gene, CDS, tRNA, and rRNA features from related GenBank records to the query assemblies using nucleotide and protein alignments, then keep the best nonredundant calls. Compare your outputs with `gateVER/*.annotation.tsv`, `gateVER/*.draft.gff3`, `gateVER/annotation_summary.json`, and `gateVER/report.md`.

For annotation integrity checks, use the gateVER-only report package in `gateVER/annotation_integrity/`. It intentionally ignores `dudleya_organelle_qc/` and evaluates only the gateVER best nonredundant annotations, cached gateVER references, and gateVER whole-genome summaries. Start with `gateVER/annotation_integrity/report.md`, then inspect `validation_checklist.tsv` for the complete procedure map. The package includes per-CDS ORF status, expected gene presence, mitochondrial core-gene presence, IR feature duplication, tRNA/rRNA plausibility, whole-genome synteny, cross-organelle identity, repeat structure, read/graph availability, and manual IGV review targets.

## Results

The chloroplast FASTA is strongly supported as a Dudleya chloroplast genome. The contig `ptg000216l_1` covers nearly all complete Dudleya chloroplast references at about 99.3 to 99.6 percent weighted nucleotide identity.

The focused `NC_085682.1` comparison in `nc085682_comparison/` compares the candidate *Dudleya setchellii* chloroplast assembly against the public *Dudleya farinosa* chloroplast reference `NC_085682.1` from https://www.ncbi.nlm.nih.gov/nuccore/NC_085682.1. It reproduces the raw BLAST result at 99.359 percent weighted identity. The downloaded NCBI files are saved in that folder as `NC_085682.1.fetched.fa` and `NC_085682.1.fetched.gb`. After removing the candidate terminal duplicate and rotating to the reference origin, the normalized projection gives 99.556017 percent similarity excluding gaps, 99.365529 percent similarity counting reference deletions, and 99.193527 percent similarity when unmapped reference bases are counted as differences. The same package reports CDS synonymous/nonsynonymous substitutions and gene-content/order comparisons.

The mitochondrial FASTA is supported as a Crassulaceae plant mitochondrial genome. The contig `ptg000317l_1` covers large portions of related mitochondrial references at about 97.3 to 97.7 percent weighted nucleotide identity. The cross-organelle comparisons found no evidence that the two FASTA labels are swapped.

The gateVER annotation-integrity pass classifies the current best nonredundant CDS calls as follows: chloroplast `51 PASS`, `25 WARN`, `9 REVIEW`; mitochondria `11 PASS`, `14 WARN`, `9 REVIEW`. It also validates expected gene content, mitochondrial core genes, chloroplast IR feature placement, noncoding RNA labels, synteny, repeat structure, and cross-organelle identity. `PASS` means the current continuous interval has a clean ATG-to-stop ORF with no internal stop. `WARN` and `REVIEW` calls remain useful homology evidence, but need boundary, split-gene, RNA-editing, or overlap review before curated annotation use.

## Visualization

Open IGV Desktop and load `igv/dudleya_organelles.igv.json` as a genome file. The chromosome dropdown should show `chloroplast` and `mitochondria`. More detailed loading notes are in `igv/README.md`.

## Limitations

These are FASTA-only analyses. Raw reads were not mapped back in this repository, and `gateVER/annotation_integrity/read_graph_tool_availability.tsv` found no read alignment, raw read, or assembly graph files under the repository root. That means coverage uniformity, SNP/indel hotspot review, and read-backed mitochondrial repeat/junction validation cannot be performed from the current repository contents. The draft annotations are homology-transfer annotations and should not be treated as curated GenBank submission annotations. The chloroplast assembly appears to include a terminal duplicate and should be de-duplicated and rotated before final submission-quality annotation.

## Citation

If you use this repository, cite the repository URL and cite the NCBI accessions used as reference evidence. The key accessions and searches are listed in `gateVER/reference_search_summary.json`, `gateVER/report.md`, and `dudleya_organelle_qc/report.md`.

For citation questions, data-use questions, or reuse requests, contact Justen Whittall at Santa Clara University (SCU) at `jwhittall@scu.edu`, or contact the repository owner at `neilsumanth9@gmail.com`.

## License

No license has been declared for this repository. Contact the repository owner before reusing or redistributing the data beyond normal GitHub viewing and citation.
