# Phylogenetic Trees

This step builds cpDNA and mtDNA phylogenetic trees from the
full callable-site consensus FASTA alignments using IQ-TREE maximum-likelihood
inference. Bootstrap support is included when
requested for the run.

## Run

- Run label: `primary`
- Method: iqtree_ml_ufboot1000

## Results

### cpDNA

- Track: `cpdna_population_sites`
- Samples: 278
- Alignment sites: 124538
- Missing bases: 347258
- Model: `GTR+F+G4`
- Tree: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.treefile`
- IQ-TREE report: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.iqtree`
- Log: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/19_bootstrap_phylogenetic_tree/cpDNA.primary.iqtree_ml.log`

### mtDNA

- Track: `mtdna_high_confidence_unique`
- Samples: 278
- Alignment sites: 44930
- Missing bases: 152444
- Model: `GTR+F+G4`
- Tree: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.treefile`
- IQ-TREE report: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.iqtree`
- Log: `/home/neil/Downloads/dudleya/Dudleya_S_genomes/full_pipeline_run/results/19_bootstrap_phylogenetic_tree/mtDNA.primary.iqtree_ml.log`
