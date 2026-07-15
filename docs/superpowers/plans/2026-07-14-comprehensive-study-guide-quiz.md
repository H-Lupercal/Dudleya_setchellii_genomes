# Comprehensive Study-Guide Quiz Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a 75-question, progressively difficult, systems-level short-answer quiz based on the existing Dudleya organelle pipeline study guide.

**Architecture:** Add one self-contained Chapter 24 Markdown file without modifying Chapters 1–23. Use three difficulty bands and continuous numbering, with early questions establishing pipeline-wide concepts and later questions requiring multi-stage debugging and scientific synthesis.

**Tech Stack:** Markdown, repository study-guide chapters, ripgrep, Perl-based numbering validation, Git.

---

### Task 1: Build the quiz content

**Files:**
- Read: `dudleya_organelle_alignment_pipeline/study_guide/01-data-flow-map.md` through `dudleya_organelle_alignment_pipeline/study_guide/23-capstone-sample-trace.md`
- Create: `dudleya_organelle_alignment_pipeline/study_guide/24-comprehensive-short-answer-quiz.md`

- [ ] **Step 1: Establish the file structure**

Create the Chapter 24 heading, a short purpose paragraph, directions asking for concise explanatory answers, and these three top-level sections:

```markdown
# Chapter 24 — Comprehensive Short-Answer Quiz

## Part I — Foundations and Pipeline Relationships (Questions 1–25)

## Part II — Application, Comparison, and Prediction (Questions 26–50)

## Part III — Debugging, Interpretation, and Synthesis (Questions 51–75)
```

Use continuous bold numbering in the existing exercise style: `**1.**` through `**75.**`.

- [ ] **Step 2: Write Questions 1–25**

Distribute the foundation questions as follows:

- 1–8: end-to-end provenance, stage dependencies, cpDNA/mtDNA separation, and initial-versus-final runs;
- 9–17: organelle haploidy/linkage, repeats, gene trees, file formats, and coordinate systems;
- 18–25: sample identity, read pairing, QC, mapping/base quality, depth versus breadth, masks, and population tracks.

Require explanations and relationships rather than one-word definitions. Do not ask for elementary Python syntax.

- [ ] **Step 3: Write Questions 26–50**

Distribute the application questions as follows:

- 26–36: trace reads through mapping, masks, haploid calling, filtering, SNP alignments, and callable consensus;
- 37–43: compare SNP-only and callable-site alignments and predict effects of changed thresholds or missing data;
- 44–50: compare phylogenetic trees, PCA, ADMIXTURE-style clustering, and Fst, including what each input and output means.

At least half of these questions must provide a small scenario, parameter change, or pair of outputs to interpret.

- [ ] **Step 4: Write Questions 51–75**

Distribute the advanced questions as follows:

- 51–60: diagnose multi-stage failures involving sample order, filenames, tracks, missingness, tool availability, or output inconsistencies;
- 61–69: evaluate scientific claims involving bootstrap support, PCA separation, K selection, Fst, linked organelle loci, and cpDNA/mtDNA discordance;
- 70–75: synthesize evidence across the full pipeline, identify the most important limitations, and formulate defensible conclusions or follow-up checks.

Ensure advanced prompts require reasoning across multiple chapters. Explicitly test the distinction between evidence, inference, and unsupported causal claims.

### Task 2: Verify and record the quiz

**Files:**
- Verify: `dudleya_organelle_alignment_pipeline/study_guide/24-comprehensive-short-answer-quiz.md`

- [ ] **Step 1: Validate numbering and question count**

Run:

```bash
perl -ne 'if (/^\*\*(\d+)\.\*\*/) { $n++; die "expected $n, found $1\n" if $1 != $n } END { die "expected 75 questions, found $n\n" unless $n == 75; print "75 questions in sequence\n" }' dudleya_organelle_alignment_pipeline/study_guide/24-comprehensive-short-answer-quiz.md
```

Expected: `75 questions in sequence`.

- [ ] **Step 2: Check format and prohibited content**

Run:

```bash
rg -n '^\s*[A-D][.)]|^\s*[1-4][.)]|Answer Key|Solutions' dudleya_organelle_alignment_pipeline/study_guide/24-comprehensive-short-answer-quiz.md
```

Expected: no output. Manually confirm every prompt is short answer, no solutions or answer-revealing hints appear, and basic Python syntax is not a substantial topic.

- [ ] **Step 3: Check topical coverage**

Run:

```bash
rg -n 'cpDNA|mtDNA|FASTQ|FASTA|BAM|BED|VCF|QC|mapping quality|depth|breadth|mask|haploid|callable|bootstrap|PCA|ADMIXTURE|Fst|discordance|tool|test' dudleya_organelle_alignment_pipeline/study_guide/24-comprehensive-short-answer-quiz.md
```

Expected: matches across all three parts. Manually compare the prompts with Chapters 1–23 and confirm the quiz does not assume untaught Stage 21 material.

- [ ] **Step 4: Check repository scope and Markdown cleanliness**

Run:

```bash
git diff --check
git status --short
```

Expected: no whitespace errors; only the new quiz file and the plan/design workflow files are involved, with no modifications to Chapters 1–23.

- [ ] **Step 5: Commit the quiz**

```bash
git add dudleya_organelle_alignment_pipeline/study_guide/24-comprehensive-short-answer-quiz.md
git commit -m "docs: add comprehensive study guide quiz"
```
