# Chapter 20 — Glossary

> Part 4 of 4 · Practice and Reference · Prev:
> [Solutions](./19-solutions.md) · Next: [Module and Function
> Index](./21-module-and-function-index.md)

Terms are alphabetical. Each entry points to the chapter where it is used in
context. Biology, bioinformatics, and software terms are mixed, because reading
this pipeline requires all three at once.

**ADMIXTURE** — a program that fits `K` ancestry clusters to genotype data and
reports a cross-validation error per K. Here it is fed pseudo-diploid, strongly
linked organelle SNPs. Its output is an exploratory model projection under
violated assumptions, not validated nuclear admixture, haplotype assignments, or
biological group counts.
[Chapters 12](./12-pca-and-clustering.md), [16](./16-pca-clustering-fst-interpretation.md).

**Base quality** — a per-base confidence score in a read. `samtools depth -q` and
`bcftools mpileup -Q` set base-quality minimums. Not the same as mapping quality.
[Chapter 4](./04-shell-and-external-tools.md).

**BED** — a region file format that is **0-based, half-open**: `start 3` means the
4th base, `end` is exclusive. The only 0-based format in the pipeline.
[Chapter 5, §5.5](./05-bioinformatics-file-formats.md).

**Breadth (of coverage)** — the fraction of a region covered by at least *N*
reads (`breadth_ge_1x`, `_5x`, `_10x`). Distinct from depth.
[Chapter 8](./08-pilot-mapping-and-investigations.md).

**Bootstrap / UFBoot** — a resampling measure of branch support in a tree.
IQ-TREE's ultrafast bootstrap (UFBoot) with `--bnni` uses ≥95 as the "strong
support" threshold, not comparable to standard bootstrap's 70.
[Chapter 15](./15-phylogenetics-interpretation.md).

**Callability** — whether a base had enough confident coverage to be trusted.
Below the minimum depth, the callable consensus writes `N`.
[Chapter 6, §6.5](./06-organelle-biology.md), [Chapter 10](./10-variants-to-alignments.md).

**Callable-site consensus** — the Stage 11 alignment covering every trusted
position: reference by default, called allele where filtered, `N` where failed or
low-depth. cpDNA 124,538 bp, mtDNA 44,930 bp.
[Chapter 10, §10.6](./10-variants-to-alignments.md).

**Comprehension / generator** — Python one-liners that build a list/dict/set
(`[... for ...]`) or lazily yield items (`(... for ...)`).
[Chapter 2, §2.6](./02-python-essentials.md).

**cpDNA (chloroplast DNA)** — the plastid genome; here ~150,274 bp with a large
inverted repeat. Treated here as one largely linked cytoplasmic locus; the
inheritance direction in *Dudleya* is not established by this repository.
[Chapter 6](./06-organelle-biology.md).

**Cross-validation (CV) error** — ADMIXTURE's measure of how well a K value
predicts held-out data; lower is better. K is chosen by lowest *mean* CV over
replicates. [Chapter 12](./12-pca-and-clustering.md).

**Dataclass** — a Python class defined by listing typed fields; `@dataclass`
generates the constructor. `frozen=True` makes it immutable.
[Chapter 2, §2.3](./02-python-essentials.md).

**Dependency injection** — passing risky, environment-dependent behavior into a
function as a parameter (with a real default) so tests can substitute a stub.
Exemplified by `check_tool(resolver=..., runner=...)`.
[Chapter 14, §14.4](./14-tool-audit.md).

**Depth** — the number of reads covering a position; `mean_depth` averages it over
a region. [Chapter 8](./08-pilot-mapping-and-investigations.md).

**Discordance** — disagreement between the cpDNA and mtDNA trees. Check resolution
and missingness before invoking biology.
[Chapter 15, §15.5](./15-phylogenetics-interpretation.md).

**FASTA** — named sequences (`>header` then bases). Used for the reference and the
alignments. [Chapter 5, §5.2](./05-bioinformatics-file-formats.md).

**FASTQ** — raw reads, four lines each. The pipeline only counts them.
[Chapter 5, §5.1](./05-bioinformatics-file-formats.md).

**Fst** — a measure of population differentiation. Here a Nei-style,
gene-diversity estimator over informative SNP sites, computed in pure Python.
[Chapter 13](./13-population-fst.md).

**Gene tree vs species tree** — a gene tree is one locus's genealogy; a species
tree is the organism's history. Organelle trees are gene trees.
[Chapter 17, §17.1](./17-uncertainty-bias-and-limits.md).

**Haploid / `--ploidy 1`** — one allele per site. Organelles are effectively
haploid, so variant calling uses `--ploidy 1`.
[Chapter 6, §6.1](./06-organelle-biology.md).

**Haplotype diversity** — the probability two random samples are different whole
sequence strings, sample-size corrected. This implementation includes `N` in the
strings, so differing missing-data patterns can inflate it.
[Chapter 13, §13.5](./13-population-fst.md).

**Heteroplasmy** — the presence of more than one organelle haplotype in an
individual; assumed away by the haploid model.
[Chapter 17, §17.2](./17-uncertainty-bias-and-limits.md).

**Inverted repeat (IR)** — two near-identical, opposite-orientation copies of a
region in the chloroplast (~`82091–107826` and `124539–150274` here). One copy is
masked for population genetics. [Chapter 6, §6.3](./06-organelle-biology.md).

**IQ-TREE** — the maximum-likelihood tree program used with `GTR+F+G4`.
[Chapter 11](./11-phylogenetic-trees.md).

**Linkage** — non-independence of sites that travel together. Each organelle is
treated as one largely linked locus, so its SNPs are not independent observations;
plant mtDNA can nevertheless recombine among repeats.
[Chapter 6, §6.2](./06-organelle-biology.md).

**Manifest** — the Stage 00 sample table; `analysis_samples.tsv` is its
authoritative primary-analysis subset.
[Chapter 7](./07-manifest-and-reference-preflight.md).

**Mapping quality (MAPQ)** — confidence that a read is placed correctly. Repeats
lower it. Permissive-MAPQ coverage overstates trustworthy mtDNA coverage.
[Chapter 5, §5.4](./05-bioinformatics-file-formats.md), [Chapter 8](./08-pilot-mapping-and-investigations.md).

**Missingness (`N` / `F_MISSING`)** — absence of a callable base. Sites with >20%
missing genotypes are filtered; PCA imputes, Fst excludes, ADMIXTURE codes `0 0`.
[Chapter 6, §6.5](./06-organelle-biology.md).

**mtDNA (mitochondrial DNA)** — the mitochondrial genome; ~243,359 bp, repeat-rich,
so only ~44,930 bp of uniquely mappable sequence is trusted.
[Chapter 6, §6.4](./06-organelle-biology.md).

**Newick** — the parenthetical tree format IQ-TREE writes; rendered with
Biopython. [Chapter 5, §5.8](./05-bioinformatics-file-formats.md).

**Nucleotide diversity** — average per-site difference across sample pairs,
counting only positions where both have real bases.
[Chapter 13, §13.5](./13-population-fst.md).

**Paralog** — a duplicated, similar sequence that can attract misplaced reads (a
mtDNA hazard). [Chapter 6, §6.4](./06-organelle-biology.md).

**PCA (principal component analysis)** — an unsupervised projection of the SNP
matrix onto axes of greatest variance; missing states are mean-imputed.
[Chapters 12](./12-pca-and-clustering.md), [16](./16-pca-clustering-fst-interpretation.md).

**PED/MAP** — PLINK text genotype format; the pipeline writes haploid calls as
homozygous "pseudo-diploid" pairs for ADMIXTURE.
[Chapter 5, §5.7](./05-bioinformatics-file-formats.md).

**Pseudo-diploid encoding** — duplicating each haploid base into a homozygous pair
so a diploid tool (ADMIXTURE) will accept it. A tooling trick, not biology.
[Chapter 5, §5.7](./05-bioinformatics-file-formats.md).

**Private variant** — an allele found in one population and no other; a set
difference in the code. [Chapter 13, §13.5](./13-population-fst.md).

**Property (`@property`)** — a method accessed like a field, recomputed on access
(e.g. `metrics.mean_depth`). [Chapter 2, §2.3](./02-python-essentials.md).

**Reference bias** — under-detection of variation in regions that differ a lot
from the single mapping reference. [Chapter 6, §6.6](./06-organelle-biology.md).

**Resumability** — reusing existing outputs unless `--force`/`--refresh-qc` is
passed, so a crashed run can restart. [Chapter 3, §3.6](./03-reusable-code-patterns.md).

**Run label** — a filename prefix (`primary`, `smoke`) threaded through every
output by `labeled_output_name`. [Chapter 3, §3.4](./03-reusable-code-patterns.md).

**SNP alignment** — the Stage 10 FASTA of variable sites only (cpDNA 2,015, mtDNA
146 columns); input to PCA, clustering, Fst.
[Chapter 10, §10.5](./10-variants-to-alignments.md).

**Subprocess / `Popen` pipe** — running external tools from Python, wiring one
tool's stdout into the next's stdin. [Chapter 3, §3.5](./03-reusable-code-patterns.md).

**Track (analysis track)** — a named BED region with a stated purpose (QC vs
population genetics). Track IDs select some downstream behavior, but the
`purpose` field itself is documentation rather than a fully enforced contract.
[Chapter 9](./09-masks-alignment-and-sample-qc.md).

**TSV** — tab-separated table with a header; the pipeline's inter-stage contract
format, read with `csv.DictReader`.
[Chapter 5, §5.9](./05-bioinformatics-file-formats.md).

**VCF** — variant call format; haploid genotypes here are single alleles (`0`,
`1`, `.`). [Chapter 5, §5.6](./05-bioinformatics-file-formats.md).

> Next: [Chapter 21 — Module and Function Index](./21-module-and-function-index.md)
