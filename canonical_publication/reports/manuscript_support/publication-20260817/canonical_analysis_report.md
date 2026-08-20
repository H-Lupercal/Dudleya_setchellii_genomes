# Canonical Dudleya organelle analysis — publication-20260817

## Status

This report is generated exclusively from the canonical dependency chain. Pre-remediation files are quarantined and were not analysis inputs.

## Sample sets

- Chloroplast: 276 QC-eligible samples.
- Mitochondria: 271 QC-eligible samples.
- Shared intersection: 271 samples, used only for concatenation.

- Incomplete read pairs excluded before analysis: 2.
- Provider-manifest entries declared absent from the deposit: 4.

- Provider manifests with an unauthenticatable self-checksum entry: 1; these manifest files remain independently SHA-256 inventoried.

- Immutable source files passing the SHA-256 inventory: 577/577.

## Read processing and evidence filters

Paired reads were adapter-trimmed and filtered with fastp, using Q20 as the qualified-base threshold, rejecting reads with more than 40% unqualified bases, and requiring length ≥50. BWA-MEM mapping included read groups. Evidence with MAPQ <20 was excluded; unmapped, secondary, supplementary, QC-failed, and duplicate records were removed. Duplicate marking used paired-read fixmate metadata, and downstream depth/pileup evidence required base quality ≥20.

Across 278 complete pairs, fastp retained 4,162,003,714/4,299,491,762 reads (96.80%) and 600,678,656,432/644,923,764,300 bases (93.14%). The base-weighted Q20 fraction increased from 98.01% to 99.18%; 643,419,232 reads (14.97%) were adapter-trimmed. The median fastp duplication-rate diagnostic was 7.90%; alignment duplicates were subsequently removed before inference.

Eligibility is organelle-specific and requires at least 80% breadth at DP≥5 over the regenerated organelle unique-mappability mask; full-reference breadth is reported separately.

## References and callable masks

- External FASTA/GenBank sequence-version consistency: 2/2 accession pairs PASS.
- Chloroplast: 150,274 bp after a self-BLAST-validated 26,690 bp redundant terminal-copy trim; 99.94% of selected query bases align to NC_085682.1 at 99.276% position-assigned HSP identity. Both IR copies are excluded only from the mappability denominator; one duplicate IR copy is excluded from population analysis.
- Mitochondria: the 243,359 bp candidate is retained intact, but only 71.11% of query bases and 42.31% of PV256627.1 are covered by qualifying alignments. Self-repeats mask 196,351 bp (80.68%); 43,182 read-supported unique bases remain in the final high-confidence mask. Among callable consensus bases, the median eligible-sample identity to the selected mitochondrial reference is 99.9883%; this mapping-conditioned concordance is not an independent assembly validation. The median repeat-mask/unique-site depth ratio is 0.005; repeat coordinates remain excluded regardless of depth. This structural discordance limits whole-mitogenome interpretation.
- The dominant external-reference HSP has reverse_complement orientation for chloroplast and same orientation for mitochondria. After assigning overlapping query positions to the highest-bitscore HSP, the same-orientation fractions among covered query bases are 54.60% and 46.00%, respectively. These are local alignment diagnostics, not evidence of global collinearity.
- External GenBank feature projection recovered 242/263 chloroplast features and 55/97 mitochondrial features. Of the projected mitochondrial features, 14/55 overlap the final read-backed high-confidence mask and 8 are fully contained. These annotations are explicitly draft projections, not de novo annotations.
- The median eligible-sample boundary/interior depth ratio is 1.432 for chloroplast and 0.000 for mitochondria; these values are evidence diagnostics, not proof of circularity.

## Variants and population statistics

Genotypes are haploid and are masked when DP <5, GQ <20, or either field is missing. Accepted biallelic SNP sites require QUAL ≥30 and no more than 20% missing genotypes. Fixed-alternate accepted sites are retained for consensus; segregating primary summaries include MAC≥1, while PCA and supplementary ADMIXTURE alone require MAC≥2. The per-input pileup depth cap is 250 reads per site.

- Chloroplast high-confidence variant sites used for consensus (including fixed alternate): 2273.
- Chloroplast primary variants (including singletons): 2261.
- Mitochondrial high-confidence variant sites used for consensus (including fixed alternate): 157.
- Mitochondrial primary variants (including singletons): 146.
- Chloroplast callable-site π range: 0–0.0029134344.
- Mitochondrial callable-site π range: 0–0.00045170257.
- Pairwise differentiation is signed Hudson ratio-of-sums FST with 1 kb block-bootstrap intervals.

Supplementary ADMIXTURE tested K=1–12 with ten fixed seeds per K. The minimum mean cross-validation error selected K=12 for chloroplast (boundary optimum) and K=12 for mitochondria (boundary optimum). These are sensitivity results under the limitations stated below.

## Phylogenetic interpretation

Separate unrooted chloroplast and mitochondrial ModelFinder trees with 1,000 SH-aLRT and 1,000 ultrafast-bootstrap replicates are primary. The partitioned concatenated tree is supplementary. The coordinate-padded chloroplast partition spans 38.2% of concatenated coordinates, but masked all-N coordinates are not treated as evidence: chloroplast contributes 71.5% of sites with at least two callable shared samples, 85.2% of jointly callable shared-sample sites, 93.9% of variable sites, including singletons, and 94.3% of parsimony-informative sites; mitochondria contributes 28.5%, 14.8%, 6.1% and 5.7%, respectively. 0 pairs of strongly supported incompatible organelle splits were detected.

## Limitations

No nuclear decoy assembly is available. MAPQ, base-quality, depth, and duplicate filters reduce but cannot eliminate NUMT/NUPT ambiguity. PCA describes linked organelle haplotype variation rather than independent loci. Diversity, differentiation, and haplotype estimates are conditional on the organelle-specific cohort and callable masks; mitochondrial estimates do not represent the repeat-rich whole candidate assembly. ADMIXTURE uses pseudo-diploid linked organelle markers and is descriptive supplementary clustering, not ancestry inference. Trees remain unrooted because no defensible outgroup was supplied.
