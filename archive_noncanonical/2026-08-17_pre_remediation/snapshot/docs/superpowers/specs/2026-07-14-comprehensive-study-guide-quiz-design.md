# Comprehensive Study-Guide Quiz Design

## Goal

Create one Markdown quiz in the existing study-guide directory that tests the
reader's understanding of the Dudleya organelle pipeline at a systems level.
The quiz should emphasize biology, bioinformatics reasoning, pipeline data flow,
method selection, interpretation, and debugging rather than elementary Python
syntax.

## Deliverable

Create:

`dudleya_organelle_alignment_pipeline/study_guide/24-comprehensive-short-answer-quiz.md`

The file will contain 75 short-answer questions and no answer key. The user will
answer questions independently and use this chat for grading or explanation.

## Structure

The quiz will progress in difficulty:

1. **Questions 1–25: Foundations and relationships.** Test core biological
   concepts, provenance, file roles, QC distinctions, and connections between
   major pipeline stages. These questions may involve recall, but should require
   explanation rather than one-word definitions.
2. **Questions 26–50: Application and comparison.** Ask the learner to trace data
   through stages, predict consequences of parameter changes, compare analytical
   methods, and explain tool or format choices.
3. **Questions 51–75: Debugging, interpretation, and synthesis.** Present
   multi-stage scenarios involving conflicting evidence, methodological limits,
   unexpected outputs, bias, and defensible biological conclusions.

Questions will be grouped under descriptive topical headings while retaining one
continuous numbering sequence.

## Coverage

The 75 questions will collectively cover:

- the end-to-end provenance chain and the separation of cpDNA and mtDNA;
- organelle inheritance, linkage, haploidy, repeats, and the gene-tree caveat;
- FASTQ, FASTA, SAM/BAM, BED, VCF, TSV, and Newick roles and coordinate rules;
- manifests, sample identity, reference validation, read pairing, and QC;
- mapping quality, base quality, depth, breadth, masks, and population tracks;
- haploid variant calling, filtering, missingness, SNP-only alignments, and
  callable-site consensus alignments;
- phylogenetic inference, branch support, bootstrap interpretation, and
  cpDNA/mtDNA discordance;
- PCA, ADMIXTURE-style clustering, pseudo-diploid encoding, and model limits;
- haplotypes, population summaries, pairwise Fst, and informative sites;
- external tools, auditability, command provenance, tests, validation, and
  dependency injection at an architectural rather than syntax-trivia level;
- uncertainty, bias, linked-locus limitations, and appropriately cautious
  scientific claims.

The quiz is based on the contents of the existing 23 study-guide chapters. It
will not silently test newly added Stage 21 material that the current guide does
not yet teach.

## Question Style

- Every item is short answer; there are no multiple-choice, true/false, or
  matching questions.
- Prompts ask for explanations, comparisons, predictions, diagnoses, or concise
  evidence-based conclusions.
- Basic questions such as defining a Python list, reproducing simple syntax, or
  naming a single obvious function are excluded.
- Some questions provide small hypothetical results or failure messages, but no
  question requires running the pipeline.
- Difficult questions may combine several chapters, but each prompt will state
  enough context to be answerable from the study guide.
- The quiz will not contain solutions, hints that disclose answers, or an answer
  key.

## Presentation

The file will begin with brief instructions explaining the intended answer
length and suggesting that the learner cite stages, files, tools, or limitations
where relevant. Section headings will make the progressive difficulty visible.
Questions will use bold continuous identifiers such as `**1.**` through
`**75.**`, following the existing exercise style.

## Quality Checks

Before completion:

- confirm exactly 75 uniquely numbered questions;
- confirm no multiple-choice options or answer key appears;
- check that elementary Python syntax is not a major focus;
- check broad coverage against Chapters 1–23;
- check that biological claims and numeric examples agree with the study guide;
- check Markdown formatting and links;
- verify that no existing study-guide chapter is modified.
