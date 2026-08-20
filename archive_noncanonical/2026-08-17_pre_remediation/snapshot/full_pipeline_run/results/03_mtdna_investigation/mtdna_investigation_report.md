# mtDNA Pilot Coverage Investigation

## Question

The first Step 3 pilot summary appeared to show low mtDNA breadth across all
pilot samples. This report checks whether that meant the mtDNA reference lacked
read support, or whether the signal was caused by mapping/filtering behavior.

## Main Finding

The original low mtDNA breadth was caused by a `samtools depth` flag
interpretation error in the Step 3 pipeline.

For `samtools depth`:

- `-q` is minimum base quality.
- `-Q` is minimum mapping quality.

The first Step 3 run used `samtools depth -q 0 -Q 20`, which reported breadth
after requiring mapping quality `>=20`. The intended permissive organelle
coverage summary is `samtools depth -q 20 -Q 0`.

After refreshing QC with the corrected command, the pilot mtDNA breadth is high:

```text
Pilot samples summarized: 15
Total cpDNA+mtDNA mapped read records: 46584669
Median input organelle mapping fraction: 0.144205
Median chloroplast breadth >=1x: 0.999993
Median mitochondrial breadth >=1x: 0.960445
Median mitochondrial breadth >=10x: 0.944646
```

Only `ABAB_MAD_LP_222_Du-589` remains flagged, and that sample is tiny
(`5378` input read records).

## What The Low-MAPQ Signal Means

The high-MAPQ-only check is still biologically useful. With MAPQ `>=20`, most
mtDNA positions lose support, while two repeat-associated blocks remain covered
across nearly all usable pilot samples:

```text
59539-84470
223117-243114
```

Those two blocks total about `44930` bp. This does not mean the rest of the
mtDNA has no reads; it means many reads outside those uniquely placed blocks
map ambiguously under the current combined cpDNA/mtDNA reference.

This matches the previous reference-verification warning that the mitochondrial
assembly has large repeats:

```text
69706 bp direct repeat: 153662-223367 vs 84115-153820
28078 bp inverted repeat: 31770-59847 vs 1-28078
```

## Chloroplast-Like Sequence Check

The reproducible mtDNA coverage is not explained by chloroplast contamination.
Local BLAST of the mtDNA reference against the normalized chloroplast reference
found only:

```text
mtDNA bp with chloroplast hit: 4651 / 243359
fraction of mtDNA with chloroplast hit: 0.019112
chloroplast-hit bp overlapping high-MAPQ consensus blocks: 873
```

## Interpretation

The mtDNA reference has broad pilot read support when mapping quality is not
used as a uniqueness filter. The scientific risk is not absence of mtDNA
coverage; it is repeat-driven ambiguous placement and possibly alternative
plant mitochondrial conformations.

For population analyses, do not treat all whole-mtDNA sites equally yet.

Recommended next mtDNA path:

1. Keep permissive MAPQ coverage metrics for sample-level mtDNA presence and
   breadth QC.
2. Add a separate MAPQ>=20 or repeat-mask track for high-confidence uniquely
   placed mtDNA variants.
3. Consider a conserved-gene mtDNA alignment for tree/PCA/Fst if whole-mtDNA
   repeat ambiguity remains high.
4. Do not use low-MAPQ/repeat regions for final SNP-based population statistics
   until a repeat/unique-placement mask is defined.
