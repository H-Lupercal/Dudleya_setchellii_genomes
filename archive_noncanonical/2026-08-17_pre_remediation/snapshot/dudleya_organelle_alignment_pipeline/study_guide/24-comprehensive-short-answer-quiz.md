# Chapter 24 — Comprehensive Short-Answer Quiz

> Part 4 of 4 · Practice and Reference · Prev: [Capstone Sample
> Trace](./23-capstone-sample-trace.md)

This quiz tests whether you can reason across the entire Dudleya organelle
pipeline, not merely recall isolated terms. Answer in complete thoughts. Most
responses should take two to five sentences; scenario and synthesis questions
may need a short paragraph. Name relevant stages, files, tools, assumptions, or
limitations whenever they strengthen your answer.

## Part I — Foundations and Pipeline Relationships (Questions 1–25)

### Pipeline architecture and provenance

**1.** Describe the major transformations that take the data from downloaded
FASTQ reads to the final phylogenetic, PCA, ADMIXTURE-style, and Fst results.
Explain why this is a provenance chain rather than a collection of unrelated
analyses.

**2.** Why must sample identity and read pairing be established before mapping,
and why must mapping and sample QC occur before haploid variant calling?

**3.** Explain why cpDNA and mtDNA are analyzed separately even though reads are
mapped to one combined organelle reference and stored in one BAM per sample.

**4.** Which upstream products are shared by several downstream analyses, and
how can a mistake in one shared product propagate into trees, PCA, clustering,
and Fst?

**5.** Why does the Stage 13 tool audit sit beside the biological data-flow chain
rather than inside it? What kind of reproducibility question does it answer?

**6.** Explain the relationship between Stages 12 and 19, Stages 14 and 20, and
Stages 16 and 18. Why should the final versions supersede the initial versions
for interpretation without requiring separate analysis modules?

**7.** What problem does a run label such as `primary` or `smoke` solve? Describe
one failure that could result from combining files produced under different run
labels.

**8.** Explain how stage summaries, reports, exact command tables, and tests
serve different but complementary roles in making the pipeline auditable.

### Organelle biology and data formats

**9.** What does it mean to treat an organelle genome as haploid, and how does
that assumption affect variant calling and genotype representation?

**10.** Why are hundreds or thousands of organelle SNPs not equivalent to the
same number of independent nuclear markers? Explain the role of physical
linkage.

**11.** Why is a cpDNA or mtDNA phylogeny a gene tree rather than a species tree?
What prevents a well-supported organelle tree from automatically establishing
the history of the species?

**12.** What is the chloroplast inverted repeat, and why does retaining both
near-identical copies in a population-genetic track risk misleading mapping and
variant calls?

**13.** Why does the pipeline distinguish permissive mtDNA coverage from the
high-confidence unique mtDNA population track? What biological signal can
repetitive sequence falsely create?

**14.** Contrast FASTQ and FASTA files in terms of what information they contain
and where each enters the pipeline.

**15.** Explain the distinct roles of a BAM file and a VCF file. What evidence is
present in the BAM that has already been summarized or discarded in the VCF?

**16.** BED is 0-based and half-open, whereas reference annotations and VCF
positions are generally 1-based. Explain why confusing these systems produces a
biological error rather than merely a cosmetic formatting error.

**17.** What information do Newick and TSV files represent in this repository?
Why is neither format interchangeable with a sequence alignment?

### Sample QC, mapping evidence, and analysis tracks

**18.** Explain how the sample manifest connects filenames, biological sample
identities, read mates, sequencing batches, and population metadata.

**19.** Why are samples missing either R1 or R2 excluded from the primary
paired-end analysis? Why would silently treating the remaining mate as an
ordinary paired sample be unsafe?

**20.** Compare the purposes of the pilot alignment and the all-sample
alignment. What questions should be settled during the pilot before scaling to
hundreds of samples?

**21.** Distinguish mapping quality from base quality. Give one example of a
poor-quality observation that each filter is meant to remove.

**22.** A region has very high mean depth but low breadth at 1×. Explain how that
pattern can occur and why mean depth alone would give a misleading picture of
callability.

**23.** Why can a broad coverage track be suitable for measuring whether a
sample contains organelle reads but unsuitable for population-genetic variant
calling?

**24.** Explain how an analysis mask converts a biological concern about repeats
or ambiguous placement into a concrete restriction on downstream computation.

**25.** The pipeline retains 275 of 280 biological sample rows for downstream
analysis. Explain why that number reflects at least two distinct kinds of QC
decisions and should not be described simply as “five samples failed mapping.”

## Part II — Application, Comparison, and Prediction (Questions 26–50)

### Mapping, masks, variants, and alignments

**26.** Trace one sample's paired FASTQ files through BWA and samtools to a
sorted, indexed BAM. For each tool, state the transformation it performs and
why the next tool needs its output.

**27.** A sample has two R1 files and two R2 files from separate lanes. Explain
why the manifest marks this differently from a simple complete pair and what
must be resolved before the standard alignment reader can use it.

**28.** A sample has strong cpDNA mapping but almost no mtDNA mapping. Give two
plausible explanations, and describe which mapping or coverage summaries you
would inspect before deciding whether to exclude it.

**29.** A trusted annotation covers bases 10 through 20 in 1-based inclusive
coordinates. Write the corresponding BED start and end coordinates, then
explain the conversion.

**30.** Suppose Stage 08 accidentally uses `cpdna_full_coverage` rather than
`cpdna_population_sites`. Predict the likely effect on the SNP set and explain
why a successful command exit would not demonstrate biological correctness.

**31.** Explain how `bcftools mpileup` and `bcftools call --ploidy 1` divide the
work of variant calling. What would be conceptually wrong about interpreting a
heterozygous diploid genotype from this model?

**32.** A position appears in the raw VCF but fails Stage 09 filtering. How is
that position represented in the Stage 11 callable consensus, and why is that
safer than automatically writing the reference base?

**33.** At one biallelic SNP the reference is `T`, the alternate allele is `C`,
and three haploid genotype calls are `1`, `0`, and `.`. State the three FASTA
bases Stage 10 should write and explain each conversion.

**34.** Predict how lowering the minimum minor-allele count and allowing a
higher missing fraction would affect the number of Stage 09 variants and the
width of the Stage 10 alignment. What new quality risks would accompany the
larger matrix?

**35.** What information is retained in a Stage 10 SNP-only alignment, and what
information about invariant or uncallable positions is absent?

**36.** Describe how Stage 11 combines the reference, filtered variants, raw-only
failed sites, the analysis track, and per-sample depth to construct a callable
consensus sequence.

### Comparing alignment types and predicting downstream effects

**37.** Why is the SNP-only alignment appropriate for PCA and Fst, while the
full callable-site alignment is the intended input for the final IQ-TREE
phylogeny?

**38.** If the minimum depth for callability increases from 1× to 5×, predict
the direction of change in Stage 11 missing bases and phylogenetic information.
Why might the stricter threshold improve reliability while reducing resolution?

**39.** If Stage 09 allows substantially more missing genotypes, explain how
that choice could affect Stage 10, PCA mean imputation, haplotype counts, and
pairwise Fst informative-site counts.

**40.** A sample is missing at every retained SNP. Explain why downstream code
should reject or exclude it rather than representing it as an ordinary point in
PCA or an ordinary individual in ADMIXTURE.

**41.** In a PCA site with observed bases `A, A, T, N`, explain how the bases are
numerically encoded and how the missing value is imputed. What assumption is
introduced by that imputation?

**42.** Why are `G,G,G,G` and `A,A,N,N` both uninformative for a biallelic PCA
matrix even though the second column contains two different characters?

**43.** Two populations have different patterns of missing calls at otherwise
similar sites. Explain how missingness could imitate biological separation in
whole-sequence haplotypes or distort the effective evidence in PCA and Fst.

### Comparing downstream analytical methods

**44.** Why does the final IQ-TREE analysis use the callable-site consensus and
a nucleotide substitution model? What does a branch length mean that a
bootstrap value does not?

**45.** Compare the purpose of the fast initial maximum-likelihood tree with the
1,000-replicate UFBoot tree using BNNI. Which question does the final analysis
answer more rigorously?

**46.** In simple terms, explain how bootstrap support is obtained by repeatedly
resampling alignment columns. What does 97% support say, and what does it not
say, about the corresponding biological grouping?

**47.** PCA places samples near one another according to variation in the input
matrix. Explain why spatial proximity on a PCA plot is not an ancestry
percentage or a direct estimate of recent gene flow.

**48.** Why does this pipeline encode haploid organelle alleles as
pseudo-diploid homozygotes for ADMIXTURE? Identify the two major model
assumptions that remain violated after this technical conversion.

**49.** Mean ADMIXTURE cross-validation error decreases continuously from K=1
through the maximum tested K=8. State the most defensible conclusion and explain
why “there are exactly eight biological populations” is not supported.

**50.** Explain how pairwise Fst is calculated conceptually from within- and
total-population diversity. Why must the `informative_sites` count accompany an
Fst estimate when comparing population pairs?

## Part III — Debugging, Interpretation, and Synthesis (Questions 51–75)

### Diagnosing multi-stage failures

**51.** Stage 11 reports that VCF sample order does not match
`included_samples.tsv`. List the most likely provenance failures and explain why
silently reordering only one file could conceal a more serious mismatch.

**52.** A rerun unexpectedly includes 276 downstream samples instead of 275.
Describe a systematic investigation using Stage 00, Stage 06, and Stage 07
outputs to identify the extra sample.

**53.** One sample has excellent `cpdna_full_coverage` breadth but poor
`cpdna_population_sites` breadth. Explain how both measurements can be correct
and what this means for sample QC versus cpDNA variant reliability.

**54.** Pilot mtDNA breadth is high when all reads are counted but collapses
under high mapping-quality requirements. Diagnose the likely genome feature and
explain why increasing sequencing depth alone may not solve it.

**55.** Stage 10 produces a FASTA whose sample names match the VCF but not the
current Stage 07 table. What should be treated as the source of truth, and which
upstream products may need to be regenerated?

**56.** Stage 09 reports 2,015 filtered cpDNA SNPs, but the Stage 10 FASTA has
2,014 columns. Give several concrete checks that can distinguish a malformed
VCF record, a stale output, and a bug in alignment construction.

**57.** IQ-TREE exits with status zero, but the expected `.treefile` is missing
or empty. Why should the pipeline still fail, and which command, prefix, log,
and report details would you inspect?

**58.** The tool audit finds `python3`, but an import check for scikit-learn
fails. Explain why the audit correctly records the dependency as missing even
though the executable exists, and why testing only `python3 --version` would be
insufficient.

**59.** A five-sample `smoke` VCF is accidentally passed to a stage using the
`primary` 275-sample metadata. Predict where validation should fail and explain
the value of failing before any biological summary is written.

**60.** A colleague changes a Stage 05 mask and reruns only Stage 17. The Fst
files are unchanged, so they conclude the mask has no effect. Explain the
provenance error and identify the downstream chain that must be regenerated.

### Auditing scientific interpretations

**61.** Two populations separate strongly on cpDNA PC1, have pairwise cpDNA Fst
of 0.45, and occupy different clades with 98 UFBoot support. Write a cautious
interpretation that combines the evidence without claiming that cpDNA proves
complete reproductive isolation.

**62.** The same two populations have mtDNA Fst of 0.00, but only three sites
were informative for that comparison. Explain why “the populations share the
same mitochondrial history” is too strong and what should be checked next.

**63.** A cpDNA branch has UFBoot support of 96, while the corresponding mtDNA
branch has support of 52. Explain what differs in evidential strength and why
the comparison should consider the very different numbers of informative
sites.

**64.** cpDNA and mtDNA trees place the same sample in conflicting groups.
Construct an ordered checklist that examines technical resolution and
missingness before invoking introgression, lineage sorting, or another
biological explanation.

**65.** cpDNA and mtDNA both support the same geographic split. Why is this
stronger than evidence from one organelle, yet still not equivalent to two
independent nuclear loci?

**66.** A report calls the organelle patterns “maternal population structure.”
Explain why the word *maternal* is not established by this repository and how
you would rewrite the claim.

**67.** ADMIXTURE selects K=8 by the lowest mean cross-validation error across
five replicates. Explain why replicate stability improves confidence in the
numerical optimization but does not repair the biological problems caused by
pseudo-diploid, linked organelle markers.

**68.** PCA clusters largely match named species. Give one defensible
interpretation and at least three reasons not to conclude that PCA has validated
the taxonomy by itself.

**69.** A high pairwise Fst is described as proof that a mountain range caused
the observed differentiation. Separate the measured result from the proposed
causal explanation and identify additional evidence needed for the causal
claim.

### Full-pipeline synthesis

**70.** For one population pair, cpDNA shows strong PCA separation, high Fst,
and a well-supported tree split, while mtDNA shows weak separation, low Fst,
and poorly supported branches. Write the strongest conclusion justified by all
of these observations together, including the one-locus and resolution caveats.

**71.** Suppose missing data is concentrated in one sequencing batch that also
corresponds mostly to one species. Explain how this confounding could influence
sample exclusion, PCA, haplotype diversity, tree placement, and Fst. Propose a
cross-stage audit to detect it.

**72.** A conservation decision requires identifying independently evolving
populations. Explain what useful evidence this organelle pipeline provides and
what key evidence it cannot supply on its own.

**73.** Explain how small fixtures, command-builder tests, output-existence
checks, and dependency injection collectively test the pipeline without
rerunning BWA, bcftools, IQ-TREE, or ADMIXTURE on the full dataset. What kinds of
biological errors can still escape these tests?

**74.** You receive only a final cpDNA tree image and are asked to reproduce it.
List the minimum provenance information and intermediate products you would
need to trace before treating the image as reproducible scientific evidence.

**75.** Summarize the pipeline as one connected argument: begin with how raw
reads become trusted organelle characters, explain how several downstream
methods examine those characters differently, name the single most important
interpretive limitation, and propose one nuclear-data follow-up that would
materially strengthen the biological conclusions.
