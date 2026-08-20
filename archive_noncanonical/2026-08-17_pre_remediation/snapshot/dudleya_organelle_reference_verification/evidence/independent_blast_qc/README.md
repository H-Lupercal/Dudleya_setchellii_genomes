# Dudleya Candidate Organelle FASTA QC and Identity Report

## Scope

This report summarizes the independent BLAST-based QC evidence generated for two raw FASTA assemblies:

- `../Dudleya_hifiasm_purged_manual_chloroplast.fa`
- `../Dudleya_hifiasm_purged_manual_mitochondria.fa`

This evidence was generated independently of the reference-verification annotation workflow. It has been retained here as a cross-check on organelle identity.

This folder contains QC and homology evidence for organelle identity. It is not a final curated annotation package and does not contain GFF3 or GTF annotation files.

## Input FASTA QC

| File | Contig | Length bp | GC percent | Ns |
|---|---:|---:|---:|---:|
| `../Dudleya_hifiasm_purged_manual_chloroplast.fa` | `ptg000216l_1` | 176,964 | 37.770 | 0 |
| `../Dudleya_hifiasm_purged_manual_mitochondria.fa` | `ptg000317l_1` | 243,359 | 45.517 | 0 |

Both FASTAs are single-contig assemblies with no ambiguous `N` bases.

## Methods

### Software

- BLAST+ 2.17.0
- `makeblastdb`
- `blastn`
- `tblastn`
- NCBI E-utilities downloads via `curl`

### Reference Sets

Chloroplast identity was tested against complete Dudleya chloroplast genomes downloaded into:

- `dudleya_chloroplast_refs.fa`

The chloroplast reference set contains:

- `NC_085682.1` Dudleya farinosa chloroplast, complete genome
- `PX244389.1` Dudleya virens chloroplast, complete genome
- `PX244390.1` Dudleya nesiotica chloroplast, complete genome
- `PX244391.1` Dudleya lanceolata chloroplast, complete genome
- `PX244392.1` Dudleya greenei chloroplast, complete genome
- `PX244393.1` Dudleya candelabrum chloroplast, complete genome
- `PX244394.1` Dudleya sp. T72600-1 chloroplast, complete genome

Mitochondrial identity was tested against complete Crassulaceae mitochondrial genomes downloaded into:

- `crassulaceae_mito_refs.fa`

The mitochondrial reference set contains:

- `PV256627.1` Graptopetalum paraguayense mitochondrion, complete genome
- `PV608516.1` Sedum sarmentosum mitochondrion, complete genome
- `PP024540.1` Rhodiola rosea mitochondrion, complete genome
- `NC_082108.1` Rhodiola juparensis mitochondrion, complete genome
- `NC_072122.1` Rhodiola tangutica mitochondrion, complete genome
- `NC_070303.1` Rhodiola crenulata mitochondrion, complete genome
- `NC_069572.1` Sedum plumbizincicola mitochondrion, complete genome

### BLAST Database Construction

The input assemblies and reference FASTAs were converted to nucleotide BLAST databases with commands equivalent to:

```bash
makeblastdb -in ../Dudleya_hifiasm_purged_manual_chloroplast.fa -dbtype nucl -out chloroplast_selfdb
makeblastdb -in ../Dudleya_hifiasm_purged_manual_mitochondria.fa -dbtype nucl -out mitochondria_selfdb
makeblastdb -in dudleya_chloroplast_refs.fa -dbtype nucl -out dudleya_chloroplast_refs_db
makeblastdb -in crassulaceae_mito_refs.fa -dbtype nucl -out crassulaceae_mito_refs_db
```

### Whole-Genome Identity Tests

Four whole-genome BLAST comparisons were performed:

| Question | Output file |
|---|---|
| Does the chloroplast-labeled FASTA match Dudleya chloroplast references? | `chloroplast_vs_cprefs.tsv` |
| Does the chloroplast-labeled FASTA instead match mitochondrial references? | `chloroplast_vs_mitorefs.tsv` |
| Does the mitochondria-labeled FASTA match Crassulaceae mitochondrial references? | `mitochondria_vs_mitorefs.tsv` |
| Does the mitochondria-labeled FASTA instead match Dudleya chloroplast references? | `mitochondria_vs_cprefs.tsv` |

The BLAST commands used `blastn`/`megablast` with tabular output containing query, reference, percent identity, alignment length, query coordinates, reference coordinates, bitscore, and e-value.

Coverage values in this report are union coverage values, so overlapping HSPs are counted once for query coverage and reference coverage.

### Marker-Gene Tests

Annotated coding sequences and proteins were downloaded from the best available complete references:

- Chloroplast marker source: `PX244389.1` Dudleya virens chloroplast
- Mitochondrial marker source: `PV256627.1` Graptopetalum paraguayense mitochondrion

Marker files:

- `PX244389_cds_na.fa`
- `PX244389_cds_aa.fa`
- `PV256627_cds_na.fa`
- `PV256627_cds_aa.fa`

Marker alignments:

- `PX244389_cds_vs_chloroplast.tsv`
- `PX244389_cds_vs_mitochondria.tsv`
- `PX244389_aa_vs_chloroplast.tsv`
- `PX244389_aa_vs_mitochondria.tsv`
- `PV256627_cds_vs_mitochondria.tsv`
- `PV256627_cds_vs_chloroplast.tsv`
- `PV256627_aa_vs_mitochondria.tsv`
- `PV256627_aa_vs_chloroplast.tsv`

For nucleotide CDS support, a strong marker hit was counted when the reference CDS had at least 80 percent query coverage and at least 70 percent nucleotide identity.

### Assembly-Shape Checks

Self-BLAST was used to identify large repeats, terminal duplicated sequence, and chloroplast-like inverted repeat structure. The self-BLAST databases are:

- `chloroplast_selfdb.*`
- `mitochondria_selfdb.*`

## Evidence: Chloroplast-Labeled FASTA

### Whole-Genome Reference Evidence

The chloroplast-labeled FASTA has near-complete high-identity matches to multiple complete Dudleya chloroplast genomes.

| Reference | Description | Query coverage percent | Reference coverage percent | Weighted identity percent | HSPs | Aligned bp |
|---|---|---:|---:|---:|---:|---:|
| `PX244389.1` | Dudleya virens chloroplast, complete genome | 100.00 | 99.99 | 99.564 | 6 | 228,749 |
| `NC_085682.1` | Dudleya farinosa chloroplast, complete genome | 99.95 | 99.83 | 99.358 | 8 | 229,162 |
| `PX244393.1` | Dudleya candelabrum chloroplast, complete genome | 100.00 | 99.97 | 99.388 | 8 | 228,845 |
| `PX244390.1` | Dudleya nesiotica chloroplast, complete genome | 100.00 | 99.99 | 99.368 | 8 | 228,867 |
| `PX244391.1` | Dudleya lanceolata chloroplast, complete genome | 99.95 | 99.99 | 99.342 | 7 | 228,752 |

Interpretation: this is strong evidence that `ptg000216l_1` is a Dudleya chloroplast genome.

### Cross-Organelle Check

The chloroplast-labeled FASTA has weaker and much less complete matches to mitochondrial references.

| Reference | Description | Query coverage percent | Weighted identity percent |
|---|---|---:|---:|
| `PV256627.1` | Graptopetalum paraguayense mitochondrion | 18.92 | 92.693 |
| `PP024540.1` | Rhodiola rosea mitochondrion | 15.35 | 93.439 |
| `NC_072122.1` | Rhodiola tangutica mitochondrion | 12.39 | 92.925 |

Interpretation: the cross-organelle signal is much weaker than the Dudleya chloroplast signal. It does not support the chloroplast-labeled FASTA being mitochondrial.

### Chloroplast Coding Marker Evidence

All 85 Dudleya virens chloroplast CDS records aligned strongly to the chloroplast-labeled FASTA. These correspond to 79 unique gene names. No Dudleya virens chloroplast CDS records aligned strongly to the mitochondria-labeled FASTA.

Selected marker genes:

| Gene | Coverage percent | Identity percent |
|---|---:|---:|
| `rbcL` | 100.0 | 99.8 |
| `matK` | 100.0 | 99.3 |
| `psbA` | 100.0 | 100.0 |
| `psaA` | 100.0 | 96.3 |
| `psbB` | 100.0 | 99.9 |
| `ndhF` | 100.0 | 99.4 |
| `ycf1` | 100.0 | 99.3 |
| `ycf2` | 100.0 | 99.9 |

Interpretation: the chloroplast-labeled FASTA contains the expected plastid protein-coding marker set.

### Chloroplast Assembly-Shape Evidence

Self-BLAST found two major repeat patterns:

| Repeat length bp | Identity percent | Query coordinates | Matching coordinates | Orientation |
|---:|---:|---|---|---|
| 26,702 | 99.944 | 150,275-176,964 | 1-26,700 | forward |
| 25,742 | 99.953 | 114,282-140,017 | 97,569-71,834 | reverse |

Interpretation:

- The 25.7 kb reverse repeat is consistent with a chloroplast inverted repeat region.
- The 26.7 kb forward terminal match suggests redundant terminal sequence in the raw assembly. If the terminal duplicate were removed, the sequence length would be approximately 150.3 kb, which matches known Dudleya plastome sizes.

## Evidence: Mitochondria-Labeled FASTA

### Whole-Genome Reference Evidence

The mitochondria-labeled FASTA has high-identity matches to multiple complete Crassulaceae mitochondrial genomes.

| Reference | Description | Query coverage percent | Reference coverage percent | Weighted identity percent | HSPs | Aligned bp |
|---|---|---:|---:|---:|---:|---:|
| `PV256627.1` | Graptopetalum paraguayense mitochondrion, complete genome | 70.89 | 42.14 | 97.697 | 107 | 174,455 |
| `NC_069572.1` | Sedum plumbizincicola mitochondrion, complete genome | 68.28 | 46.82 | 97.300 | 107 | 170,311 |
| `PV608516.1` | Sedum sarmentosum mitochondrion, complete genome | 65.82 | 61.37 | 97.431 | 114 | 163,924 |
| `NC_072122.1` | Rhodiola tangutica mitochondrion, complete genome | 63.71 | 36.35 | 97.495 | 149 | 161,424 |
| `PP024540.1` | Rhodiola rosea mitochondrion, complete genome | 55.48 | 36.11 | 97.664 | 166 | 160,388 |

Interpretation: this supports `ptg000317l_1` as a Crassulaceae plant mitochondrial assembly. The support is family-level rather than Dudleya-specific because a complete Dudleya mitochondrial genome was not available in this reference set.

### Cross-Organelle Check

The mitochondria-labeled FASTA has only tiny matches to Dudleya chloroplast references.

| Reference | Description | Query coverage percent | Weighted identity percent |
|---|---|---:|---:|
| `NC_085682.1` | Dudleya farinosa chloroplast | 1.24 | 77.676 |
| `PX244393.1` | Dudleya candelabrum chloroplast | 1.24 | 77.635 |
| `PX244390.1` | Dudleya nesiotica chloroplast | 1.24 | 77.635 |
| `PX244389.1` | Dudleya virens chloroplast | 1.24 | 77.635 |
| `PX244391.1` | Dudleya lanceolata chloroplast | 1.24 | 77.601 |

Interpretation: this strongly argues against the mitochondria-labeled FASTA being a mislabeled chloroplast genome.

### Mitochondrial Coding Marker Evidence

Of 31 Graptopetalum paraguayense mitochondrial CDS records, 23 aligned strongly to the mitochondria-labeled FASTA. No mitochondrial CDS records aligned strongly to the chloroplast-labeled FASTA under the same nucleotide threshold.

Strongly supported mitochondrial genes included:

`ATP1`, `ATP4`, `ATP6`, `ATP8`, `ATP9`, `COX1`, `COX2`, `COX3`, `CYTB`, `ND2`, `ND4`, `ND4L`, `ND5`, `ND7`, `ND9`, `ccmB`, `ccmC`, `ccmFc`, `matR`, `rpl10`, `rpl5`, `rps14`, and `rps7`.

Selected marker genes:

| Gene | Coverage percent | Identity percent |
|---|---:|---:|
| `COX1` | 100.0 | 99.2 |
| `COX2` | 100.0 | 99.0 |
| `COX3` | 100.0 | 98.7 |
| `CYTB` | 100.0 | 99.4 |
| `ATP1` | 100.0 | 99.2 |
| `ATP6` | 99.7 | 99.2 |
| `ND2` | 100.0 | 99.3 |
| `ND4` | 100.0 | 99.6 |
| `ND5` | 99.0 | 99.7 |
| `ND7` | 100.0 | 99.7 |
| `matR` | 100.0 | 99.4 |
| `ccmB` | 100.0 | 99.7 |

Interpretation: the mitochondria-labeled FASTA contains a strong set of plant mitochondrial coding markers.

### Mitochondrial Assembly-Shape Evidence

Self-BLAST found large internal repeats:

| Repeat length bp | Identity percent | Query coordinates | Matching coordinates | Orientation |
|---:|---:|---|---|---|
| 69,706 | 99.999 | 153,662-223,367 | 84,115-153,820 | forward |
| 28,078 | 100.000 | 31,770-59,847 | 28,078-1 | reverse |

Interpretation:

- Large repeats are common in plant mitochondrial assemblies.
- FASTA-only evidence cannot determine whether these repeats reflect true mitochondrial structure, alternate conformations, or assembly redundancy.
- Read-backed validation would be needed to resolve repeat structure and molecule configuration.

## Bottom-Line Calls

| FASTA | Identity call | Evidence strength | Important caveat |
|---|---|---|---|
| `../Dudleya_hifiasm_purged_manual_chloroplast.fa` | Dudleya chloroplast genome | Strong | Raw assembly appears to include a terminal duplicate; trim/rotate before final reference annotation or submission. |
| `../Dudleya_hifiasm_purged_manual_mitochondria.fa` | Crassulaceae plant mitochondrial genome | Strong at family/organelle level | No complete Dudleya mitochondrial reference was used; large repeats need read-backed validation. |

There is no evidence from these analyses that the two FASTA labels are swapped.

## Key Files in This Folder

| File | Purpose |
|---|---|
| `report.md` | This report. |
| `chloroplast_vs_cprefs.tsv` | Whole-genome BLAST of chloroplast FASTA against Dudleya chloroplast references. |
| `chloroplast_vs_mitorefs.tsv` | Cross-organelle BLAST of chloroplast FASTA against mitochondrial references. |
| `mitochondria_vs_mitorefs.tsv` | Whole-genome BLAST of mitochondrial FASTA against Crassulaceae mitochondrial references. |
| `mitochondria_vs_cprefs.tsv` | Cross-organelle BLAST of mitochondrial FASTA against Dudleya chloroplast references. |
| `dudleya_chloroplast_refs.fa` | Dudleya chloroplast reference FASTAs. |
| `crassulaceae_mito_refs.fa` | Crassulaceae mitochondrial reference FASTAs. |
| `PX244389_cds_na.fa` | Dudleya virens chloroplast nucleotide CDS markers. |
| `PX244389_cds_aa.fa` | Dudleya virens chloroplast protein markers. |
| `PV256627_cds_na.fa` | Graptopetalum paraguayense mitochondrial nucleotide CDS markers. |
| `PV256627_cds_aa.fa` | Graptopetalum paraguayense mitochondrial protein markers. |
| `PX244389_cds_vs_chloroplast.tsv` | Chloroplast CDS marker evidence for the chloroplast FASTA. |
| `PV256627_cds_vs_mitochondria.tsv` | Mitochondrial CDS marker evidence for the mitochondrial FASTA. |

The complete folder contains 50 generated QC/reference/BLAST files.

## Limitations

- This is FASTA-only evidence. No raw reads were mapped back to the assemblies.
- This report confirms organelle identity and major QC observations, but it does not produce a final curated annotation.
- Plant mitochondrial genomes can be rearranged, repeat-rich, and structurally multipartite. The mitochondrial FASTA should not be interpreted as a finished molecule without read-backed support.
- The chloroplast FASTA should be de-duplicated and rotated before final chloroplast annotation or GenBank-style submission.
