# Phylogenetic Trees

This step builds cpDNA and mtDNA phylogenetic trees from the Step 10
full callable-site consensus FASTA alignments using IQ-TREE
maximum-likelihood inference. This run uses a fast ML search and does
not yet include bootstrap support.

## Run

- Run label: `primary`
- Method: IQ-TREE maximum-likelihood fast search

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Samples: 275
- Alignment sites: 124538
- Missing bases: 127485
- Model: `GTR+F+G4`
- Tree: `dudleya_organelle_alignment_pipeline/results/12_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile`
- IQ-TREE report: `dudleya_organelle_alignment_pipeline/results/12_phylogenetic_tree/cpDNA.primary.iqtree_ml.iqtree`
- Log: `dudleya_organelle_alignment_pipeline/results/12_phylogenetic_tree/cpDNA.primary.iqtree_ml.log`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Samples: 275
- Alignment sites: 44930
- Missing bases: 31313
- Model: `GTR+F+G4`
- Tree: `dudleya_organelle_alignment_pipeline/results/12_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile`
- IQ-TREE report: `dudleya_organelle_alignment_pipeline/results/12_phylogenetic_tree/mtDNA.primary.iqtree_ml.iqtree`
- Log: `dudleya_organelle_alignment_pipeline/results/12_phylogenetic_tree/mtDNA.primary.iqtree_ml.log`
