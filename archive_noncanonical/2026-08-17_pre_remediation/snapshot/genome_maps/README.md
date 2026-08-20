# Circular Organelle Genome Maps

This folder contains circular maps of the *Dudleya setchellii* chloroplast and
mitochondrial reference assemblies.

These are genome-structure figures, not sample phylogenies. They use the
reference FASTA files and transferred GFF3 annotations from
`dudleya_organelle_reference_verification/`.

## Outputs

| File | What it shows |
|---|---|
| [cpDNA.circular_genome_map.png](cpDNA.circular_genome_map.png) | Chloroplast circular genome map with CDS strand tracks, tRNA/rRNA features, SNP-density ring, and GC-content deviation. |
| [cpDNA.circular_genome_map.svg](cpDNA.circular_genome_map.svg) | Editable vector version of the chloroplast map. |
| [cpDNA.circular_genome_map.pdf](cpDNA.circular_genome_map.pdf) | PDF version of the chloroplast map. |
| [mtDNA.circular_genome_map.png](mtDNA.circular_genome_map.png) | Mitochondrial circular genome map with CDS strand tracks, tRNA/rRNA features, SNP-density ring, and GC-content deviation. |
| [mtDNA.circular_genome_map.svg](mtDNA.circular_genome_map.svg) | Editable vector version of the mitochondrial map. |
| [mtDNA.circular_genome_map.pdf](mtDNA.circular_genome_map.pdf) | PDF version of the mitochondrial map. |
| [genome_map_summary.tsv](genome_map_summary.tsv) | Input files, sequence lengths, annotation counts, GC fraction, and output filenames. |

## Inputs

| Organelle | FASTA | GFF3 |
|---|---|---|
| cpDNA | `dudleya_organelle_reference_verification/references/chloroplast.normalized.fa` | `dudleya_organelle_reference_verification/annotations/chloroplast.gff3` |
| mtDNA | `dudleya_organelle_reference_verification/references/mitochondria.fa` | `dudleya_organelle_reference_verification/annotations/mitochondria.gff3` |

The maps were generated with [build_circular_genome_maps.py](build_circular_genome_maps.py)
using `pycirclize`, matplotlib, and Biopython in the local `.tools/bioconda-env`
environment.
