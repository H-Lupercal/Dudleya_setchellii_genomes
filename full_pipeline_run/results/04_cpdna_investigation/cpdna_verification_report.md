# cpDNA Pilot Coverage Verification

## Question

This verification checks whether the chloroplast pilot alignments have broad
read support, and whether any low mapping-quality regions reflect biological
chloroplast repeat structure or a more concerning cross-organelle/reference
problem.

## Main Finding

The cpDNA reference has broad pilot read support. Excluding the tiny
`ABAB_MAD_LP_222_Du-589` sample, permissive coverage with
`samtools depth -q 20 -Q 0` is effectively complete:

```text
Usable pilot samples: 14
Median chloroplast breadth >=1x: 0.999996
Minimum chloroplast breadth >=1x among usable samples: 0.999887
Median chloroplast breadth >=10x across all pilot samples: 0.999834
```

The lower high-MAPQ breadth is expected chloroplast inverted-repeat behavior,
not evidence that cpDNA is poorly supported.

## High-MAPQ Placement Check

With a stricter unique-placement check, `samtools depth -q 0 -Q 20`, the median
chloroplast breadth at `>=1x` across usable pilot samples is:

```text
0.753284
```

This means about one quarter of the chloroplast reference does not have
consistently unique high-MAPQ placement. Self-BLAST shows why.

## Normalized Chloroplast Repeat Structure

Self-BLAST of the normalized chloroplast mapping reference found one major
reverse repeat pair:

```text
82091-107826    vs    124539-150274
length: 25742 bp
identity: 99.953%
orientation: reverse
```

Merged large-repeat sequence covers:

```text
51472 bp
0.342521 of the chloroplast reference
```

The high-MAPQ non-consensus portion of the chloroplast reference is:

```text
46345 bp
```

The large repeat explains:

```text
46323 / 46345 bp
0.999525 of the high-MAPQ non-consensus sequence
```

Interpretation: nearly all cpDNA high-MAPQ ambiguity is the expected duplicated
inverted repeat, not a sample failure.

## Cross-Organelle Check

Local BLAST of the normalized chloroplast reference against the mitochondrial
reference found limited similarity:

```text
cpDNA bp with mtDNA hit: 4119 / 150274
fraction of cpDNA with mtDNA hit: 0.027410
mtDNA-hit bp overlapping high-MAPQ non-consensus sequence: 2654
fraction of high-MAPQ non-consensus explained by mtDNA hits: 0.057266
```

Interpretation: cross-organelle similarity exists but is not the main reason
for cpDNA low-MAPQ placement. The inverted repeat is the dominant driver.

## Annotation Implications

Feature overlap with the repeat/high-MAPQ regions is consistent with standard
chloroplast structure:

```text
CDS: 85 total; 18 mostly repeat-overlapping; 69 mostly high-MAPQ unique
gene: 119 total; 30 mostly repeat-overlapping; 91 mostly high-MAPQ unique
rRNA: 4 total; 4 mostly repeat-overlapping; 0 mostly high-MAPQ unique
tRNA: 30 total; 10 mostly repeat-overlapping; 20 mostly high-MAPQ unique
```

The rRNA annotations fall in the repeated region, as expected for chloroplast
inverted-repeat structure.

## Recommendation

The cpDNA pilot verification supports moving forward with chloroplast
all-sample processing.

For downstream population genetics:

1. Use permissive MAPQ coverage for sample-level cpDNA breadth/depth QC.
2. Use a repeat-aware SNP mask for PCA, Fst, structure/admixture-style plots,
   and tree-building.
3. Either mask one copy of the inverted repeat or collapse IR-equivalent sites
   so duplicated IR sequence is not counted twice as independent evidence.
4. Treat the normalized chloroplast reference as mapping-ready, but keep the IR
   intervals documented:

```text
IR copy 1: 82091-107826
IR copy 2: 124539-150274
```

Unlike mtDNA, cpDNA does not need a separate conserved-gene fallback path before
the next all-sample pilot/full-run step.
