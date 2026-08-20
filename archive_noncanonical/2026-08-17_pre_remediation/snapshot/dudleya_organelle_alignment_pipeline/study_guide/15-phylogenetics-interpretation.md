# Chapter 15 — Reading the Phylogenetic Trees

> Part 3 of 4 · Interpretation · Prev: [The Tool
> Audit](./14-tool-audit.md) · Next: [Reading PCA, Clustering, and
> Fst](./16-pca-clustering-fst-interpretation.md)

You know how the trees are built ([Chapter 11](./11-phylogenetic-trees.md)). This
chapter is about reading them without over-claiming. The recurring theme: a tree
from this pipeline is a *gene tree of one linked organelle locus*, and every
interpretive move has to respect that.

## 15.1 What the tree is a tree *of*

The Stage 19 tree is the maximum-likelihood genealogy of 275 organelle
haplotypes — one per sampled individual — under the `GTR+F+G4` model, with 1,000
ultrafast bootstrap (UFBoot) support values. `[RESULT]` Read that literally:

- The tips are individuals, labeled by sample ID.
- The topology is the best-supported branching order of their *organelle*
  sequences.
- Because each organelle is one largely linked cytoplasmic locus ([Chapter 6,
  §6.2](./06-organelle-biology.md)), the tree is a single locus's history, not the
  species' history. A clade means “these organelle sequences share a recent
  common ancestor,” which is not the same as “these individuals form a species.”
  The direction of organelle inheritance has not been established for *Dudleya*
  in this repository. `[BIO]`

## 15.2 Branch lengths and what "close" means

Branch lengths are in **substitutions per site** (the `--seqtype DNA`,
`GTR+F+G4` model estimates them). Two tips joined by short branches have few
inferred organelle substitutions between them; long branches mean more. Because
the alignment is the callable-site consensus, a long branch reflects inferred
substitutions at observed bases; IQ-TREE treats `N` as unknown, not as a
substitution. Concentrated missingness can weaken placement or support but does
not directly add branch length. Treat a very long isolated branch as a prompt to
check mapping quality, reference distance, contamination, and callable fraction
before believing it is biologically distinct. The [IQ-TREE
FAQ](https://iqtree.github.io/doc/Frequently-Asked-Questions) documents that gaps
and `N` are treated as unknown characters. `[BIO]`

## 15.3 Bootstrap support: UFBoot, and the ≥95 convention

Each internal branch in the Stage 19/20 tree carries a UFBoot value from 0–100.
UFBoot with the `--bnni` correction is designed so that **≥95 indicates strong
support** for that branch (the convention differs from standard bootstrap, where
70 is the rough threshold — UFBoot values are not directly comparable to
standard bootstrap percentages). `[BIO]` Practical reading:

- **≥95**: the data strongly support this branch *for this locus and model*.
- **<95**: unresolved; do not build a story on it.
- A well-supported *deep* split still only tells you about organelle history.

The Stage 12 first-pass tree has **no** support values at all — it is a topology
sketch. Any support-based claim must come from Stage 19/20. Confusing the two is
the most common misreading, and [Chapter 11, §11.7](./11-phylogenetic-trees.md)
flags it. `[BIO]`

## 15.4 The trees are unrooted

IQ-TREE here is run without a specified outgroup, so the inferred trees are
**unrooted** — the drawing may look rooted, but the placement of the "root" in
the figure is an artifact of rendering, not an inference. `[CODE]` Do not read
directionality ("A gave rise to B") or a most-ancestral lineage from an unrooted
tree. If you need a rooted interpretation, that requires an outgroup or a rooting
method applied outside this pipeline, and any such rooted figure elsewhere in the
repository should be treated as a separate analysis with its own assumptions.

## 15.5 cpDNA versus mtDNA: reading discordance

The pipeline deliberately builds trees for two separate cytoplasmic markers.
Comparing them is informative, but their evidence is not automatically
independent because organelles can share transmission history.

**When the two trees agree**, that is concordance: two organelle markers recovered
the same sequence grouping, although shared inheritance may help produce the
agreement. **When they disagree**
(discordance), the possibilities, in the order you should check them:

1. **Resolution.** mtDNA has only ~146 informative SNPs versus cpDNA's ~2,015
   ([Chapter 10](./10-variants-to-alignments.md)), so the mtDNA tree is simply
   less resolved. Much apparent discordance is the mtDNA tree being unable to
   support the split, not contradicting it — check the mtDNA UFBoot values first.
2. **Sample-specific callability.** Overall mtDNA consensus missingness is only
   0.2534% (31,313 of 275 × 44,930 cells), but a clade involving samples with
   unusually low callable fractions still deserves scrutiny. `[RESULT]`
3. **Only then, biology.** Genuine cpDNA/mtDNA discordance can reflect real
   organelle history — for example different capture/introgression histories of
   the two organelles, or incomplete lineage sorting. But this is the *last*
   explanation to reach for, never the first. `[BIO]`

## 15.6 What you may and may not say

| You may say | You may **not** say |
|---|---|
| "These individuals share a well-supported cpDNA clade (UFBoot ≥95)." | "These individuals are a species / are more closely related overall." |
| "The cpDNA and mtDNA trees agree on grouping X." | "The organelle trees establish the nuclear relationships." |
| "The mtDNA tree is poorly resolved here (low UFBoot, few SNPs)." | "The mtDNA tree contradicts the cpDNA tree" (before ruling out resolution). |
| "Branch lengths suggest sample Y is divergent for cpDNA." | "Sample Y is a distinct taxon" (from one unrooted locus). |
| "This is the final bootstrap-supported tree (Stage 19/20)." | "The Stage 12 topology is confident." |

## 15.7 A short reading checklist

Before you make any tree claim, confirm:

- [ ] You are looking at the Stage 19/20 (bootstrap) tree, not Stage 12/14.
- [ ] The branch you care about has UFBoot ≥95.
- [ ] You are not reading direction from an unrooted tree.
- [ ] For discordance, you have checked mtDNA resolution and missingness first.
- [ ] Your claim is about *organelle-haplotype* history, tagged as such.

Reading trees well is mostly discipline about what a single-locus, unrooted,
model-based estimate can carry. The interpretation of PCA, clustering, and Fst
follows the same discipline in [Chapter
16](./16-pca-clustering-fst-interpretation.md), and the full accounting of
uncertainty is [Chapter 17](./17-uncertainty-bias-and-limits.md).

> Next: [Chapter 16 — Reading PCA, Clustering, and Fst](./16-pca-clustering-fst-interpretation.md)
