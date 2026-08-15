# Phylogenetic Trees

This step builds cpDNA and mtDNA phylogenetic trees from the
full callable-site consensus FASTA alignments using IQ-TREE maximum-likelihood
inference. Bootstrap support is included when
requested for the run.

## Run

- Run label: `primary`
- Method: iqtree_ml_ufboot10000

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Samples: 275
- Alignment sites: 124538
- Missing bases: 127485
- Model: `GTR+F+G4`
- Tree: `dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree_10000/cpDNA.primary.iqtree_ml.treefile`
- IQ-TREE report: `dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree_10000/cpDNA.primary.iqtree_ml.iqtree`
- Log: `dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree_10000/cpDNA.primary.iqtree_ml.log`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Samples: 275
- Alignment sites: 44930
- Missing bases: 31313
- Model: `GTR+F+G4`
- Tree: `dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree_10000/mtDNA.primary.iqtree_ml.treefile`
- IQ-TREE report: `dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree_10000/mtDNA.primary.iqtree_ml.iqtree`
- Log: `dudleya_organelle_alignment_pipeline/results/19_bootstrap_phylogenetic_tree_10000/mtDNA.primary.iqtree_ml.log`
