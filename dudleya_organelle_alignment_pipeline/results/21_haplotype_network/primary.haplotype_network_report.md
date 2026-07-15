# Stage 21 Haplotype Networks (primary)

This stage applies complete-case filtering to the Stage 10 haploid SNP alignments, then builds separate cpDNA and mtDNA networks with `pegas::haploNet`.

Nodes are haplotypes, node area is sample frequency, colored sectors are species groups, and edge labels are mutational steps. These are descriptions of sequence relationships and haplotype sharing, not ancestry proportions.

## cpDNA

- Samples: 275
- Complete SNP sites: 1977 of 2015 (dropped 38)
- Haplotypes: 151
- Network edges: 6708
- Figure: `dudleya_organelle_alignment_pipeline/results/21_haplotype_network/cpDNA.primary.haplotype_network.png`
- Assignments: `dudleya_organelle_alignment_pipeline/results/21_haplotype_network/cpDNA.primary.haplotype_assignments.tsv`
- PopART NEXUS: `dudleya_organelle_alignment_pipeline/results/21_haplotype_network/cpDNA.primary.popart.nex`

## mtDNA

- Samples: 275
- Complete SNP sites: 116 of 146 (dropped 30)
- Haplotypes: 58
- Network edges: 803
- Figure: `dudleya_organelle_alignment_pipeline/results/21_haplotype_network/mtDNA.primary.haplotype_network.png`
- Assignments: `dudleya_organelle_alignment_pipeline/results/21_haplotype_network/mtDNA.primary.haplotype_assignments.tsv`
- PopART NEXUS: `dudleya_organelle_alignment_pipeline/results/21_haplotype_network/mtDNA.primary.popart.nex`

## Interpretation limits

Organelle sites are linked and represent a single nonrecombining lineage per organelle. The networks do not estimate nuclear population structure, admixture fractions, direction of gene flow, or the timing of shared ancestry.
