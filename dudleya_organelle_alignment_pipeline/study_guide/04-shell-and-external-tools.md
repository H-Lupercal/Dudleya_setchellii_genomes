# Chapter 4 — The Shell and External Tools

> Part 1 of 4 · Foundations · Prev: [Reusable Code
> Patterns](./03-reusable-code-patterns.md) · Next: [Bioinformatics File
> Formats](./05-bioinformatics-file-formats.md)

The Python in this pipeline is mostly *glue*: it decides which reads, regions,
and samples to use, then hands the heavy lifting to standard bioinformatics
tools. This chapter explains how that hand-off works and what each tool does, at
the depth you need to read the code. The formats these tools consume and produce
are the next chapter's job; here we focus on the commands.

## 4.1 The tool environment and the `PATH` prefix

The tools live in a local, git-ignored conda-style environment at
`.tools/bioconda-env/`. Tool-dependent commands are run with that directory
prepended to `PATH`:

```bash
env PATH="$PWD/.tools/bioconda-env/bin:$PATH" \
  python3 dudleya_organelle_alignment_pipeline/scripts/run_pilot_alignment.py
```

`env PATH="...:$PATH" <command>` runs `<command>` with a modified `PATH` for that
one invocation. Because the pipeline calls tools by bare name (`bwa`,
`samtools`), whichever `bwa` is first on `PATH` wins — so the prefix guarantees
the environment's pinned versions are used. The pipeline finds tools with
`shutil.which("bwa")`, which is Python's equivalent of the shell's `which`: it
returns the resolved path or `None`. Stages that need a tool call a
`require_*` helper first and raise a clear error if it is missing, rather than
letting a cryptic "command not found" surface deep inside a subprocess. `[CODE]`

The versions pinned in this project are listed in [`../README.md`](../README.md)
and audited by Stage 13 ([Chapter 14](./14-tool-audit.md)). The ones that matter
for reading the code: `bwa` 0.7.19, `samtools`/`bcftools` 1.23, `IQ-TREE` 3.1.2,
`PLINK` 1.9, `ADMIXTURE` 1.3.0.

## 4.2 Mapping: `bwa mem | samtools view | samtools sort`

The core alignment for one sample is a three-stage shell pipe, built in
[`pilot_alignment.py`](../pilot_alignment.py) `run_alignment_commands` and
streamed with `subprocess.Popen` ([Chapter 3, §3.5](./03-reusable-code-patterns.md)):

```text
bwa mem -t 4 dudleya_cp_mt.fa reads_R1.fastq.gz reads_R2.fastq.gz
  | samtools view -@ 4 -b -F 4 -q <min_mapq> -
  | samtools sort -@ 4 -o sample.tmp.bam -
```

Reading it left to right:

- **`bwa mem`** aligns the paired reads to the combined cpDNA+mtDNA reference and
  emits alignment records (SAM) on standard output. `-t 4` is threads.
- **`samtools view -b -F 4 -q <min_mapq> -`** reads that stream (`-` means
  stdin), outputs BAM (`-b`), and filters: `-F 4` drops unmapped reads (bit 4 in
  the SAM flag), and `-q <min_mapq>` drops reads below a mapping-quality
  threshold. In the pipeline the mapping default is `min_mapq = 0`, so at this
  step all *mapped* reads are kept. `[CODE]`
- **`samtools sort -o sample.tmp.bam -`** sorts by coordinate into a temporary
  BAM, which is then atomically renamed to the final name once the pipe
  succeeds.

Writing to a temp file and renaming only on success is a deliberate safety
move: a half-written BAM never masquerades as a finished one, which keeps the
resumability check in [Chapter 3, §3.6](./03-reusable-code-patterns.md) honest.

### Post-alignment QC commands

After the BAM exists, `run_qc_commands` runs four more `samtools` commands:

- **`samtools index`** builds the `.bai` index so regions can be fetched
  randomly.
- **`samtools flagstat`** reports mapped/paired/duplicate counts.
- **`samtools idxstats`** reports mapped-read counts *per reference record* —
  this is how the pipeline separates chloroplast reads from mitochondria reads,
  since both are records in the same BAM.
- **`samtools depth`** reports per-base coverage, which drives all the breadth
  and mean-depth metrics.

## 4.3 The `-q` versus `-Q` trap in `samtools depth`

This is the single most important command-line subtlety in the pipeline, and the
code comments it explicitly. In `samtools depth`, the two quality flags are
*backwards* from what you might guess:

```python
def build_depth_command(bam_path: Path, min_mapq: int, min_baseq: int) -> list[str]:
    # In `samtools depth`, `-q` means minimum base quality and
    # `-Q` means minimum mapping quality.
    return ["samtools", "depth", "-aa", "-q", str(min_baseq), "-Q", str(min_mapq),
            bam_path.as_posix()]
```

So lowercase `-q` is **base** quality and uppercase `-Q` is **mapping** quality —
the opposite of the mnemonic most people expect. A test locks this down exactly,
asserting the flag order for `min_mapq=7, min_baseq=19`. `[TEST]` `-aa` forces
*all* positions to be reported, including zero-coverage ones, so that breadth
denominators are the full reference length rather than only the covered
positions. Getting this wrong once, earlier in the project's history, is what the
Stage 02 report means by "after correcting the `samtools depth` quality flags"
([Chapter 8](./08-pilot-mapping-and-investigations.md)).

## 4.4 Variant calling: `bcftools mpileup | bcftools call`

Stage 08 genotypes samples with a two-command `bcftools` pipe, built in
[`variant_calling.py`](../variant_calling.py) `build_bcftools_commands`:

```text
bcftools mpileup -Ou --threads 4 --ignore-RG --max-depth 10000
    -q <min_mapq> -Q <min_baseq> -a FORMAT/DP,FORMAT/AD
    -f dudleya_cp_mt.fa -R <track>.bed -b <bam_list>.txt
  | bcftools call --threads 4 --ploidy 1 -m -v -Oz -o out.vcf.gz
```

- **`bcftools mpileup`** stacks the reads from every BAM in `<bam_list>.txt` at
  each position and summarizes the pileup. Here the flag meanings are the
  conventional ones: `-q` is minimum **mapping** quality and `-Q` is minimum
  **base** quality (note this is the reverse of `samtools depth` in §4.3 — the
  two tools genuinely differ, which is why the pipeline is careful about both).
  `-R <track>.bed` restricts calling to the population-genetic region for that
  organelle, `--max-depth 10000` caps pileup depth per sample, `--ignore-RG`
  treats each BAM as one sample regardless of read groups, and
  `-a FORMAT/DP,FORMAT/AD` adds per-sample depth and allele-depth annotations.
  `-Ou` emits uncompressed BCF to stream efficiently into the next command.
- **`bcftools call --ploidy 1 -m -v`** does the genotype calling. **`--ploidy 1`
  is the heart of the whole pipeline's biology**: organelle genomes are haploid,
  so each sample gets one allele, not a diploid genotype. `-m` is the
  multiallelic caller, and `-v` outputs *variant sites only*. `-Oz` writes
  bgzipped VCF. `[CODE]`

Two cleanup commands follow: **`bcftools reheader -N <sample_names>`** replaces
the auto-generated sample names (which would otherwise be BAM file paths) with
the real sample IDs, and **`bcftools index -t`** builds the `.tbi` index. The
defaults for calling are stricter than for mapping: `min_mapq = 20`,
`min_baseq = 20`, so variants come only from confidently placed, confidently
read bases. `[CODE]`

## 4.5 Variant filtering: `bcftools view`

Stage 09 filters the raw calls down to clean biallelic SNPs, built in
[`variant_filtering.py`](../variant_filtering.py):

```text
bcftools view --threads 4 -m2 -M2 -v snps
    --min-ac 2:minor -i 'F_MISSING<=0.2' -Oz -o filtered.vcf.gz raw.vcf.gz
```

- **`-m2 -M2`** keeps sites with a minimum and maximum of 2 alleles — i.e.
  strictly biallelic.
- **`-v snps`** keeps only SNPs (drops indels).
- **`--min-ac 2:minor`** requires the *minor* allele to appear at least twice, so
  singletons (a variant seen in exactly one sample, often an error) are dropped.
- **`-i 'F_MISSING<=0.2'`** keeps only sites genotyped in at least 80% of
  samples (`F_MISSING` is the fraction of missing genotypes). `[CODE]`

A test asserts every one of these flags is present, so the filter definition is
pinned. `[TEST]` These four choices — biallelic, SNP-only, no singletons,
≤20% missing — are the exact line between "raw calls" and "sites we trust for
population genetics," and their downstream effect (2,475 → 2,015 cpDNA;
190 → 146 mtDNA) is a `[RESULT]` you can read in the Stage 09 summary.

## 4.6 Trees: IQ-TREE

Stage 12/19 infer maximum-likelihood trees with IQ-TREE, built in
[`phylogenetic_tree.py`](../phylogenetic_tree.py):

```text
iqtree -s <alignment>.fa --seqtype DNA -m GTR+F+G4 --prefix <out>
    -T 4 --safe --redo --quiet [--fast] [-B 1000 --bnni]
```

- **`-m GTR+F+G4`** is the substitution model: GTR (general time-reversible)
  with empirical base frequencies (`+F`) and gamma-distributed rate variation
  across sites (`+G4`, four rate categories). This is a standard, flexible DNA
  model.
- **`--fast`** (Stage 12 only) runs a quick heuristic search — good enough to see
  the topology.
- **`-B 1000 --bnni`** (Stage 19 only) runs 1,000 ultrafast bootstrap replicates
  with the BNNI optimization to reduce bootstrap overestimation. This is the
  final, publication-grade run.
- `--redo` overwrites prior output, `--safe` uses safe numerical mode, `--quiet`
  suppresses console chatter.

The two mutually exclusive modes are controlled in Python:
`fast = not args.full_search and not args.bootstrap_replicates`. So requesting
bootstraps automatically turns *off* `--fast`, which is exactly what
distinguishes Stage 12 from Stage 19. Two tests assert `--fast` appears for the
fast run and `-B 1000 --bnni` (and *not* `--fast`) for the bootstrap run.
`[TEST]`

## 4.7 Clustering: PLINK + ADMIXTURE

Stage 16/18 cluster samples with ADMIXTURE, which needs PLINK binary input.
Built in [`admixture_analysis.py`](../admixture_analysis.py):

```text
plink --file <prefix> --make-bed --out <prefix>          # PED/MAP -> BED/BIM/FAM
admixture --cv --seed=<seed> -j4 <prefix>.bed <K>         # cluster into K groups
```

- **`plink --make-bed`** converts the human-readable PED/MAP text
  ([Chapter 5](./05-bioinformatics-file-formats.md)) into PLINK's binary
  `.bed/.bim/.fam` triple that ADMIXTURE reads.
- **`admixture --cv <prefix>.bed K`** fits `K` ancestry clusters and reports a
  cross-validation (CV) error. The pipeline sweeps `K = 1..8`, runs five seeded
  replicates per K in the final Stage 18, and picks the K with the lowest *mean*
  CV error. `--seed=<seed>` makes each replicate reproducible; the code advances
  the seed per replicate (`seed + replicate - 1`). `[CODE]`

Because ADMIXTURE runs in a working directory and writes files named after its
input, the pipeline runs it with `cwd=output_dir` and then renames the generic
`.Q`/`.P` outputs to include the K and replicate. This bookkeeping is the bulk of
the module and is covered in [Chapter 12](./12-pca-and-clustering.md).

## 4.8 Tools the pipeline runs *in Python*, not via subprocess

Not every "tool" is an external binary. Several stages do their statistics in
Python libraries:

- **PCA** uses `scikit-learn`'s `PCA` and `numpy`
  ([`pca_analysis.py`](../pca_analysis.py)).
- **Figures** use `matplotlib` (PCA scatter, structure bars, tree rendering).
- **Tree rendering** parses Newick with `Bio.Phylo` from Biopython
  ([`tree_visualization.py`](../tree_visualization.py)).
- **Fst and diversity** are computed with plain Python and `collections.Counter`
  — no external popgen tool at all ([`population_genetics.py`](../population_genetics.py)).

These still count as dependencies, which is why the tool audit checks for
`python_matplotlib`, `python_sklearn`, `python_biopython`, and friends by
importing them, not just by finding an executable ([Chapter 14](./14-tool-audit.md)).

### `MPLCONFIGDIR`

Matplotlib wants a writable config/cache directory. The plotting stages set
`os.environ.setdefault("MPLCONFIGDIR", "/tmp/dudleya_matplotlib")` and the run
commands export the same variable. `setdefault` means "use this only if the
variable is not already set," so your environment wins if you have configured
matplotlib elsewhere. If you ever see a matplotlib permissions warning, this is
the knob. `[CODE]`

## 4.9 Command provenance: the `commands.tsv` files

Mapping (02, 06), variant calling (08), filtering (09), and trees (12, 19) write
a `commands.tsv` recording exact invocations, produced by joining each argument
list with `shlex_join` (which quotes anything containing spaces). ADMIXTURE also
shells out, but records its commands in `primary.admixture_summary.tsv` rather
than a dedicated command table. The tool audit runs version probes and records
their results, not a complete command log. Python-only stages document key
parameters in their reports. Command provenance is therefore extensive but its
location and completeness vary by stage. `[CODE]`
`[CODE]`

> Next: [Chapter 5 — Bioinformatics File Formats](./05-bioinformatics-file-formats.md)
