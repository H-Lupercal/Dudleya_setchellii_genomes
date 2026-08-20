# Reference And Pilot Preflight

This step validates the combined cpDNA/mtDNA reference, records tool
availability, prepares indexes only when tools are installed, and writes
a representative pilot sample table. It does not align reads.

## Reference

- Reference: `dudleya_organelle_reference_verification/references/dudleya_cp_mt.fa`
- chloroplast: 150274 bp (PASS; expected 150274)
- mitochondria: 243359 bp (PASS; expected 243359)

## Tool Availability

- bwa: FOUND (`/home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/bwa`)
- samtools: FOUND (`/home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/samtools`)
- fastp: FOUND (`/home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/fastp`)
- fastqc: FOUND (`/home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/fastqc`)
- multiqc: FOUND (`/home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/multiqc`)
- bcftools: FOUND (`/home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/bcftools`)

## Index Status

- samtools_faidx: CREATED_OR_UPDATED. Created with samtools faidx.
- bwa: CREATED_OR_UPDATED. Created with bwa index.

## Pilot Sample Set

- Pilot samples selected: 15
- Source table: `dudleya_organelle_alignment_pipeline/results/00_manifest/analysis_samples.tsv`
- Missing-mate samples are not eligible for this pilot set.
