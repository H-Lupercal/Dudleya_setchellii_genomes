# Bioinformatics Tool Audit

This report records whether the local tools needed for the Dudleya
cpDNA/mtDNA workflow are installed and visible on `PATH`.

## Summary

- Audit label: `primary`
- Tools checked: 25
- Tools found: 23
- Missing required current-pipeline tools: none
- Missing required remaining-goal tools: none
- Missing recommended tools: snakemake, python_ete3

## Tool Checks

| Tool | Necessity | Status | Version | Path | Use |
|---|---|---|---|---|---|
| python3 | required_current | FOUND | Python 3.14.6 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/python3 | pipeline scripts, tests, PCA/Fst helper code |
| bwa | required_current | FOUND | bwa | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/bwa | read mapping to combined cpDNA/mtDNA reference |
| samtools | required_current | FOUND | samtools 1.23.1 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/samtools | BAM sorting/indexing/depth/QC |
| bcftools | required_current | FOUND | bcftools 1.24 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/bcftools | haploid variant calling, filtering, and VCF indexing |
| fastp | required_current | FOUND | fastp 1.3.6 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/fastp | read QC/trimming if rerunning raw-read QC |
| fastqc | required_current | FOUND | FastQC v0.12.1 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/fastqc | read QC reports |
| multiqc | required_current | FOUND | multiqc, version 1.35 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/multiqc | aggregate QC reports |
| iqtree | required_current | FOUND | IQ-TREE version 3.1.2 for Linux x86 64-bit built May 17 2026 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/iqtree | maximum-likelihood cpDNA/mtDNA phylogenetic trees |
| FastTree | recommended_remaining | FOUND | FastTree 2.2.0 Double precision: | /home/neil/miniconda3/bin/FastTree | quick approximate ML tree checks |
| plink | required_remaining | FOUND | PLINK v1.9.0-b.8 64-bit (22 Oct 2024) | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/plink | PCA matrix handling and ADMIXTURE input preparation |
| admixture | required_remaining | FOUND | ****                   ADMIXTURE Version 1.3.0                  **** | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/admixture | structure/admixture-style clustering and empirical K selection |
| vcftools | recommended_remaining | FOUND | VCFtools (0.1.17) | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/vcftools | VCF-based population summary checks |
| bedtools | recommended_remaining | FOUND | bedtools v2.31.1 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/bedtools | interval QC and mask cross-checking |
| Rscript | required_remaining | FOUND | Rscript (R) version 4.3.3 (2024-02-29) | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/Rscript | PCA/tree/admixture/Fst plotting outputs |
| snakemake | recommended_remaining | MISSING |  |  | optional integration with Snakemake orchestration |
| python_matplotlib | required_remaining | FOUND | 3.10.9 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/python3 | PCA scatterplots, tree rendering, and static PNG/PDF figures |
| python_pandas | required_remaining | FOUND | 3.0.3 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/python3 | figure-ready metadata tables and plotting data frames |
| python_sklearn | required_remaining | FOUND | 1.9.0 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/python3 | PCA calculation and variance summaries |
| python_biopython | required_remaining | FOUND | 1.87 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/python3 | Newick tree parsing and scripted tree rendering |
| python_seaborn | recommended_remaining | FOUND | 0.13.2 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/python3 | polished statistical plots |
| python_ete3 | recommended_remaining | MISSING |  | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/python3 | alternate Newick tree visualization |
| r_ggplot2 | required_remaining | FOUND | 3.5.2 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/Rscript | R-based PCA/admixture/Fst plots |
| r_ape | required_remaining | FOUND | 5.8.1 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/Rscript | R-based phylogenetic tree parsing/rendering |
| r_pegas | required_current | FOUND | 1.3 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/Rscript | haploid cpDNA/mtDNA haplotype networks |
| r_patchwork | recommended_remaining | FOUND | 1.3.2 | /home/neil/Downloads/dudleya/Dudleya_S_genomes/.tools/bioconda-env/bin/Rscript | combining multiple R figure panels |

## Interpretation

- The completed mapping, variant, consensus, and first-pass tree steps can be reproduced with the required current-pipeline tools.
- The remaining planned analyses have their required external tools available.
