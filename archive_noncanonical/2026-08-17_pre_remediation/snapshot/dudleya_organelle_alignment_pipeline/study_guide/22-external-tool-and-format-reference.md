# Chapter 22 — External Tool and File-Format Reference

> Part 4 of 4 · Practice and Reference · Prev:
> [Module and Function Index](./21-module-and-function-index.md) · Next:
> [Capstone](./23-capstone-sample-trace.md)

A lookup reference for the external tools and file formats. The teaching versions
are [Chapter 4](./04-shell-and-external-tools.md) (tools) and [Chapter
5](./05-bioinformatics-file-formats.md) (formats); this chapter is the quick
table you keep open while reading code.

## 22.1 External tools

Versions are the versions observed by Stage 13 and documented in
[`../README.md`](../README.md). The environment specification uses a mixture of
exact versions and version ranges, so these are not all strictly pinned.

| Tool | Version | Used for | Where in code |
|---|---|---|---|
| `bwa` | 0.7.19 | map reads (`bwa mem`) | `pilot_alignment.run_alignment_commands` |
| `samtools` | 1.23 | view/sort/index/flagstat/idxstats/depth | `pilot_alignment` |
| `bcftools` | 1.23 | mpileup, call, view (filter), reheader, index | `variant_calling`, `variant_filtering` |
| `iqtree`/`iqtree2` | 3.1.2 | ML trees | `phylogenetic_tree.build_iqtree_command` |
| `plink` | 1.9 | PED/MAP → BED/BIM/FAM | `admixture_analysis.run_plink_make_bed` |
| `admixture` | 1.3.0 | K clustering + CV | `admixture_analysis.run_admixture_for_k` |
| `fastp`, `fastqc`, `multiqc` | 1.3.5 / 0.12.1 / 1.35 | read QC (audited; not driven in these stages) | `tool_audit`, Stage 01 checks |
| `vcftools`, `bedtools` | 0.1.17 / 2.31.1 | recommended cross-checks | `tool_audit` only |
| Python `numpy`, `scikit-learn` | — | PCA matrix + decomposition | `pca_analysis` |
| Python `matplotlib` | 3.11.0 | all figures | `pca_analysis`, `admixture_analysis`, `tree_visualization` |
| Python `Bio` (Biopython) | 1.87 | Newick parsing/drawing | `tree_visualization` |
| `Rscript` + `ggplot2`/`ape`/`patchwork` | 4.5.3 / … | audited R plotting stack | `tool_audit` only |

Run tool-dependent stages with the environment on `PATH`:
`env PATH="$PWD/.tools/bioconda-env/bin:$PATH" python3 …`
([Chapter 4, §4.1](./04-shell-and-external-tools.md)).

## 22.2 Exact commands the pipeline runs

Copy-read these against the code; each is built as an argument list, never a
shell string.

**Mapping (Stage 02/06)** — `pilot_alignment.run_alignment_commands`:

```text
bwa mem -t <threads> <ref> <r1> <r2>
  | samtools view -@ <threads> -b -F 4 -q <min_mapq> -
  | samtools sort -@ <threads> -o <sample>.tmp.bam -
samtools index <bam>
samtools flagstat <bam>          # -> qc/<sample>.flagstat.txt
samtools idxstats <bam>          # -> qc/<sample>.idxstats.tsv
samtools depth -aa -q <min_baseq> -Q <min_mapq> <bam>   # -> qc/<sample>.depth.tsv
```

**Variant calling (Stage 08)** — `variant_calling.build_bcftools_commands`:

```text
bcftools mpileup -Ou --threads <t> --ignore-RG --max-depth 10000
    -q <min_mapq=20> -Q <min_baseq=20> -a FORMAT/DP,FORMAT/AD
    -f <ref> -R <track>.bed -b <bam_list>.txt
  | bcftools call --threads <t> --ploidy 1 -m -v -Oz -o <pre_reheader>.vcf.gz
bcftools reheader -N <sample_names>.txt -o <final>.vcf.gz <pre_reheader>.vcf.gz
bcftools index -t <final>.vcf.gz
```

**Variant filtering (Stage 09)** — `variant_filtering.build_filter_command`:

```text
bcftools view --threads <t> -m2 -M2 -v snps
    --min-ac 2:minor -i 'F_MISSING<=0.2' -Oz -o <filtered>.vcf.gz <raw>.vcf.gz
bcftools index -t <filtered>.vcf.gz
```

**Trees (Stage 12/19)** — `phylogenetic_tree.build_iqtree_command`:

```text
iqtree -s <aln>.fa --seqtype DNA -m GTR+F+G4 --prefix <out>
    -T <t> --safe --redo --quiet [--fast] [-B 1000 --bnni]
```

**Clustering (Stage 16/18)** — `admixture_analysis`:

```text
plink --file <prefix> --make-bed --out <prefix>
admixture --cv --seed=<seed> -j<t> <prefix>.bed <K>
```

## 22.3 The two quality-flag conventions

The single most confusable detail. The two tools disagree, and the pipeline is
careful about both ([Chapter 4, §4.3](./04-shell-and-external-tools.md)):

| Flag | `samtools depth` | `bcftools mpileup` |
|---|---|---|
| `-q` | minimum **base** quality | minimum **mapping** quality |
| `-Q` | minimum **mapping** quality | minimum **base** quality |

## 22.4 File formats

| Format | Extension | Coordinates | Read/written by |
|---|---|---|---|
| FASTQ | `.fastq.gz` | — | Stage 00 (names), 02/06 (counts) |
| FASTA (reference) | `.fa` | 1-based | Stage 01, 02, 06, 08, 11 |
| FASTA index | `.fa.fai` | 1-based lengths | Stage 02, 06 |
| SAM/BAM | `.bam` | 1-based | Stage 02/06 write, 08 reads |
| BED | `.bed` | **0-based half-open** | Stage 05 writes; 06, 08, 11 read |
| VCF | `.vcf.gz` | 1-based `POS` | Stage 08/09 write; 09, 10, 11 read |
| FASTA alignment | `.fa` | shared column index | Stage 10/11 write; 12/19, 15, 16/18, 17 read |
| PLINK PED/MAP, BED/BIM/FAM | `.ped/.map`, `.bed/.bim/.fam` | marker index | Stage 16/18 |
| Newick | `.treefile` | branch lengths | Stage 12/19 write; 14/20 read |
| TSV | `.tsv` | column-keyed | every stage |

## 22.5 Key genotype and coordinate conventions

- **Haploid genotype** in VCF: `0` = REF, `1` = ALT, `.` = missing (no `0/1`).
  `genotype_to_base` decodes to a base or `N`.
- **BED ↔ 1-based**: BED `start` = 1-based start − 1; BED `end` = 1-based
  inclusive end. `interval_to_bed_fields` and its inverse.
- **Pseudo-diploid PED**: each haploid base duplicated (`A → A A`), missing → `0
  0`, with six leading PLINK columns (`fid iid 0 0 0 -9`).
- **SAM FLAG bit 4** = unmapped; `samtools view -F 4` drops it.
- **Run label** prefixes outputs: `primary.` (deliverable) or `smoke.` (test),
  via `labeled_output_name`.

## 22.6 Output filename patterns

Once you know the run label and organelle, you can predict any output name:

| Pattern | Example |
|---|---|
| `<organelle>.<label>.raw.vcf.gz` | `cpDNA.primary.raw.vcf.gz` |
| `<organelle>.<label>.filtered.vcf.gz` | `mtDNA.primary.filtered.vcf.gz` |
| `<organelle>.<label>.snp_alignment.fa` | `cpDNA.primary.snp_alignment.fa` |
| `<organelle>.<label>.callable_consensus.fa` | `mtDNA.primary.callable_consensus.fa` |
| `<organelle>.<label>.iqtree_ml.treefile` | `cpDNA.primary.iqtree_ml.treefile` |
| `<organelle>.<label>.pca.png` | `mtDNA.primary.pca.png` |
| `<organelle>.<label>.bestK<K>.structure.png` | `cpDNA.primary.bestK8.structure.png` |
| `<label>.<stage>_summary.tsv` / `_report.md` | `primary.variant_filtering_summary.tsv` |

> Next: [Chapter 23 — Capstone: Tracing One Sample End to End](./23-capstone-sample-trace.md)
