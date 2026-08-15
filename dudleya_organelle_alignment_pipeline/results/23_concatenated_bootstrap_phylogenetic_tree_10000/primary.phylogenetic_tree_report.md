# Phylogenetic Trees

This step builds cpDNA and mtDNA phylogenetic trees from the
full callable-site consensus FASTA alignments using IQ-TREE maximum-likelihood
inference. Bootstrap support is included when
requested for the run.

## Run

- Run label: `primary`
- Method: iqtree_ml_ufboot10000

## Results

### cpDNA_mtDNA

- Track: `cpdna_then_mtdna`
- Samples: 275
- Alignment sites: 169468
- Missing bases: 158798
- Model: `GTR+F+G4`
- Tree: `dudleya_organelle_alignment_pipeline/results/23_concatenated_bootstrap_phylogenetic_tree_10000/cpDNA_mtDNA.primary.iqtree_ml.treefile`
- IQ-TREE report: `dudleya_organelle_alignment_pipeline/results/23_concatenated_bootstrap_phylogenetic_tree_10000/cpDNA_mtDNA.primary.iqtree_ml.iqtree`
- Log: `dudleya_organelle_alignment_pipeline/results/23_concatenated_bootstrap_phylogenetic_tree_10000/cpDNA_mtDNA.primary.iqtree_ml.log`
