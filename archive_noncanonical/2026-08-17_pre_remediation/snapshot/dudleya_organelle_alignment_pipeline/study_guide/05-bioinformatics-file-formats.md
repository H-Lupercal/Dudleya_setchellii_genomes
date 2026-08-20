# Chapter 5 — Bioinformatics File Formats

> Part 1 of 4 · Foundations · Prev: [Shell and
> External Tools](./04-shell-and-external-tools.md) · Next: [Organelle
> Biology](./06-organelle-biology.md)

Every arrow in the data-flow map is a file in one of about nine formats. You do
not need format-standard expertise; you need enough to read the pipeline's
parsers and reason about what a file contains. This chapter gives you exactly
that, and it dwells on the one thing the pipeline is most careful about:
coordinate systems.

## 5.1 FASTQ — raw reads

FASTQ holds sequencing reads. Each read is four lines: an `@` header, the
bases, a `+` separator, and per-base quality characters. The pipeline never
parses read *content*; it only *counts* reads, by counting lines and dividing by
four, in [`pilot_alignment.py`](../pilot_alignment.py):

```python
def count_fastq_records(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    ...
    if line_count % 4 != 0:
        raise AlignmentError(f"FASTQ line count is not divisible by 4: {path}")
    return line_count // 4
```

Two lessons here. First, FASTQ files are usually gzipped (`.fastq.gz`), so the
code picks `gzip.open` or `open` by suffix — a pattern you will see again for
VCFs. Second, a line count not divisible by four means a truncated or corrupt
file, so the code raises rather than reporting a wrong count. `[CODE]` The read
count becomes the denominator of the "what fraction of input reads mapped to an
organelle" metric.

The *filenames* carry all the sample metadata, and parsing them is Stage 00's
entire job. An Illumina name like
`CY_RED_LP_202_Du-561_S192_L005_R1_001.fastq.gz` encodes population code, plant
IDs, sequencing sample, lane, and read direction. That parsing is
[Chapter 7](./07-manifest-and-reference-preflight.md).

## 5.2 FASTA — reference and alignments

FASTA holds named sequences: a `>` header line, then the sequence. The combined
organelle reference is two records in one file:

```text
>chloroplast
ACGT...            (150,274 bases across many wrapped lines)
>mitochondria
ACGT...            (243,359 bases)
```

The pipeline reads FASTA with hand-written parsers (no Biopython needed for
this) in several modules; the pattern is always "on a `>` line, start a new
record; otherwise append bases." From
[`prepare_reference_and_pilot.py`](../prepare_reference_and_pilot.py):

```python
def read_fasta_lengths(fasta_path: Path) -> dict[str, int]:
    ...
    if line.startswith(">"):
        current_name = line[1:].split()[0]   # record name = first token after ">"
    else:
        current_length += len(line)
```

`line[1:].split()[0]` takes everything after `>` and keeps the first
whitespace-delimited token as the record name — so `>chloroplast description`
becomes `chloroplast`. The pipeline uses FASTA three ways: the input reference,
the SNP-only alignments (Stage 10), and the full callable-site consensus
alignments (Stage 11). In an *alignment* FASTA, every record is the same length
and the position index is shared across samples — that is what makes it an
alignment rather than just a bag of sequences. `[CODE]`

## 5.3 `.fai` — the FASTA index

`samtools faidx reference.fa` writes `reference.fa.fai`, a small text table with
one row per record: name, length, and byte offsets. The pipeline reads only the
first two columns to learn each organelle's length without re-scanning the whole
FASTA:

```python
def read_fai_lengths(fai_path: Path) -> dict[str, int]:
    ...
    lengths[fields[0]] = int(fields[1])
```

So `chloroplast → 150274`, `mitochondria → 243359`. The mapping stage requires
this file (and the `bwa index` files) to exist before it will run, via
`require_reference_indexes` — missing indexes raise a clear "re-run reference
preflight" error rather than a confusing tool crash. `[CODE]`

## 5.4 SAM/BAM — aligned reads

SAM is the text format for aligned reads; BAM is its compressed binary form. The
pipeline works with BAM and never parses it directly in Python — it lets
`samtools` summarize it (§4.2). You still need two SAM concepts to read the
code:

- **The FLAG** is a bitfield describing each read. Bit 4 (value `4`) means
  "unmapped." That is why `samtools view -F 4` drops unmapped reads: `-F` means
  "exclude reads with this bit set." `[CODE]`
- **MAPQ**, mapping quality, is a per-read confidence that the read is placed
  correctly. Higher is better. It is central to organelle work because repeats
  create reads that *could* map to two places; those get low MAPQ. The
  distinction between "permissive MAPQ" and "high MAPQ" coverage is the whole
  reason the mtDNA region had to be investigated ([Chapter
  8](./08-pilot-mapping-and-investigations.md)).

The per-BAM summaries the pipeline actually consumes are `idxstats` (mapped reads
per record) and `depth` (coverage per position), both plain TSV.

## 5.5 BED — regions, and the coordinate trap

BED files list genomic intervals: record, start, end, optional name. This is the
format the analysis masks (Stage 05) are written in, and it hides the single
most error-prone idea in genome coordinates.

**BED is 0-based and half-open.** A BED interval `chloroplast 0 3` covers
positions 1, 2, 3 in ordinary 1-based counting — it starts *at* index 0 and
stops *before* index 3. Most biologists think in **1-based inclusive**
coordinates, where `1–3` means positions 1, 2, 3. The pipeline converts between
them in exactly one place, and a test pins the conversion:

```python
def interval_to_bed_fields(record, start_1based, end_1based, name):
    # 1-based inclusive -> 0-based half-open
    return [record, str(start_1based - 1), str(end_1based), name]
```

So the 1-based inclusive interval `82091–107826` becomes the BED fields
`82090 107826`. `[CODE]` A unit test asserts exactly this
(`interval_to_bed_fields(..., 82091, 107826, ...) == ["chloroplast", "82090",
"107826", ...]`). `[TEST]` To keep everyone honest, Stage 05 also writes
`analysis_regions.tsv` recording *both* coordinate systems side by side, so you
never have to guess which one a number is in.

When code reads a BED interval back, it reverses the conversion: BED
`start_0based` plus 1 gives the 1-based start, and `end_0based` is already the
1-based inclusive end. You will see `start_1based=bed_start + 1,
end_1based=bed_end` in [`all_sample_alignment.py`](../all_sample_alignment.py)
and `range(start_0based + 1, end_0based + 1)` in
[`callable_consensus.py`](../callable_consensus.py). Whenever you compare a BED
number to a position, check which base you are in — this is where off-by-one bugs
breed, and the pipeline's discipline about it is not decoration.

## 5.6 VCF — variants and genotypes

VCF (Variant Call Format) holds called variants. After `##` header lines, one
`#CHROM` column line names the samples, then each data row is one variant site.
The columns before the samples are `CHROM POS ID REF ALT QUAL FILTER INFO
FORMAT`, and everything from column 10 on is one genotype per sample.

For this haploid pipeline, a genotype is a single allele index: `0` = the REF
allele, `1` = the ALT allele, `.` = missing. (Diploid VCFs would show `0/1`;
here `--ploidy 1` produces single values.) The genotype may carry extra
fields after a colon, like `0:12` (allele 0, depth 12). The pipeline decodes a
genotype to a base with a function that appears, identically, in both
[`snp_alignment.py`](../snp_alignment.py) and
[`callable_consensus.py`](../callable_consensus.py):

```python
def genotype_to_base(genotype_field: str, ref: str, alt: str) -> str:
    genotype = genotype_field.split(":", 1)[0]        # drop the ":12" part
    allele = genotype.replace("|", "/").split("/", 1)[0]  # first allele
    if allele == "0":
        return ref.upper()
    if allele == "1":
        return alt.upper()
    return "N"                                        # missing -> N
```

`split(":", 1)[0]` keeps only the genotype before any extra fields;
`replace("|", "/").split("/", 1)[0]` takes the first allele whether the
separator is `|` or `/`. So `0:12 → REF`, `1:9 → ALT`, `.:0 → N`. A test feeds a
tiny three-sample VCF and asserts the decoded sequences are exactly `"AC"`,
`"GT"`, and `"NN"`. `[TEST]` The pipeline reads gzipped VCFs (`.vcf.gz`) with the
same suffix-based `gzip.open` trick as FASTQ.

## 5.7 PED/MAP — PLINK genotype text

ADMIXTURE cannot read a FASTA, so Stage 16/18 rewrites the SNP alignment as a
PLINK PED/MAP pair. The MAP file lists markers (one per SNP); the PED file has
one row per sample with two alleles per marker. Because organelle DNA is
haploid, the pipeline duplicates each called base into a homozygous pair — this
is the "pseudo-diploid" encoding. From
[`admixture_analysis.py`](../admixture_analysis.py):

```python
for base in sequence:
    allele = base if base in BASES else "0"   # A/C/G/T kept; N -> "0" (missing)
    genotype_fields.extend([allele, allele])   # duplicate: "A" -> "A A"
ped_handle.write(" ".join([sample_id, sample_id, "0", "0", "0", "-9",
                           *genotype_fields]) + "\n")
```

The six leading fields are PLINK's family/individual/parents/sex/phenotype
columns (`-9` = missing phenotype). A test asserts a sample with sequence `TG`
becomes `... -9 T T G G` and a missing base becomes `0 0`. `[TEST]` Keep the word
"pseudo-diploid" in mind: it is a *tooling* trick to satisfy a diploid program,
not a biological claim. Together with strong linkage among organelle SNPs, it
means the resulting plots are exploratory outputs under violated model
assumptions—not validated admixture or haplotype-assignment estimates ([Chapter
16](./16-pca-clustering-fst-interpretation.md)).

## 5.8 Newick — trees

IQ-TREE writes trees in Newick format: nested parentheses with tip labels and
branch lengths, ending in a semicolon:

```text
(DU-1:0.1,(DU-2:0.2,DU-3:0.3):0.4);
```

Each `name:length` is a tip and its branch length; parentheses group clades;
the numbers after a closing paren are internal branch lengths (and, in the
bootstrap trees, support values). The pipeline never parses Newick by hand — it
hands the `.treefile` to Biopython's `Bio.Phylo.read(path, "newick")` for
rendering ([`tree_visualization.py`](../tree_visualization.py)). What you read in
the tree is [Chapter 15](./15-phylogenetics-interpretation.md).

## 5.9 TSV — the connective tissue

Every table the pipeline writes — manifests, summaries, QC decisions, PCA
coordinates, Fst matrices — is a tab-separated file with a header row. TSV is
the pipeline's lingua franca precisely because it is trivial to read with
`csv.DictReader`, diff in git, and eyeball in a terminal. The summary TSVs are
the contracts between stages ([Chapter 3, §3.3](./03-reusable-code-patterns.md)):
each stage's `read_*_inputs` opens the previous stage's summary TSV and trusts a
specific set of columns. If you want to know precisely what a stage consumes,
that column list is the answer.

## 5.10 Format-to-stage cheat sheet

| Format | Written by | Read by | Coordinate system |
|---|---|---|---|
| FASTQ (`.fastq.gz`) | sequencer | Stage 00 (names), 02/06 (mapping, counts) | n/a |
| FASTA reference | reference project | Stage 01, 02, 06, 08, 11 | 1-based |
| `.fai` | `samtools faidx` (Stage 01) | Stage 02, 06 | 1-based lengths |
| BAM (`.bam`) | Stage 02/06 | Stage 08 (mpileup) | 1-based |
| BED (`.bed`) | Stage 05 | Stage 06, 08, 11 | **0-based half-open** |
| VCF (`.vcf.gz`) | Stage 08/09 | Stage 09, 10, 11 | 1-based `POS` |
| FASTA alignment | Stage 10/11 | Stage 12/19 (trees), 15, 16/18, 17 | shared column index |
| PED/MAP, BED/BIM/FAM | Stage 16/18 | ADMIXTURE | marker index |
| Newick (`.treefile`) | Stage 12/19 | Stage 14/20 | branch lengths |
| TSV | every stage | the next stage | column-keyed |

Carry the BED row in your head: it is the only format in the pipeline that is
0-based, and confusing it with the 1-based positions in FASTA, BAM, and VCF is
the classic organelle-coordinates mistake.

> Next: [Chapter 6 — Organelle Biology That Changes the Code](./06-organelle-biology.md)
