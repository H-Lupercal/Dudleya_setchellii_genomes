# Dudleya Conservation Genomics Pipeline

Dr. Justen Whittall, Evan Hackstadt, Karina Martinez, Dante Cable

Santa Clara University (SCU)

*Original README last updated: 6/16/2026. Local source/credit note updated:
6/28/2026.*

## Source and Credit

This directory is a local copy of the Dudleya conservation genomics pipeline
originally published at:

https://github.com/evanhackstadt/dudleya

Pipeline credit belongs to Dr. Justen Whittall, Evan Hackstadt, Karina Martinez,
and Dante Cable at Santa Clara University. This copy is included here so the
pipeline can be inspected and adapted alongside the Dudleya organelle reference
work in the parent repository.


## About

This pipeline is used to study the conservation genomics of *Dudleya setchellii*, an endangered succulent endemic to Santa Clara county, California. The research project is led by Dr. Justen Whittall (Biology Department, Santa Clara University). Preliminary results were presented at the 2026 West Coast Biological Sciences Undergraduate Research Conference. Journal publication is anticipated to be forthcoming.

The results should be reproducible, but note that the pipeline is specifically designed for use on SCU's WAVE High Performance Cluster (HPC) and our specific conda environments and datasets.


## Pipeline Overview

This pipeline takes a collection of whole-genome DNA sequences, cleans them, maps them to a reference genome, and performs population analyses. It is robust to low coverage sequencing data and preserves uncertainty for the PCA and distance matrix.

The pipeline produces three primary outputs: 
* PCA plot
* Genetic distance matrix
* Phylogenetic tree.

### Key Pipeline Steps
Steps 1-3 are performed on every sample. Then, all samples are aggregated by ANGSD for steps 4-9.
1. fastp - filtering & QC
2. BWA MEM - alignment to reference genome --> .bam
3. Samtools - sort, index, and get QC report
4. ANGSD - population analysis --> .beagle, .mafs, .bcf
5. bcftools - convert .bcf --> .vcf --> .phy
6. IQTREE - produce phylogenetic tree
7. PCAngsd - produce PCA plot
8. ngsDist - produce genetic distance matrix
9. Custom aggregate QC report


## Pipeline Architecture

`Snakefile` contains the pipeline itself. It is a list of rules (terminal commands / scripts) to process and analyze the data.

`config/` contains two types of **config files** needed to run the pipeline:
* The **profile** (`config.yaml`)
  * Contains parameters telling Snakemake to use SLURM and specifying default resources
  * MUST be named config.yaml and shouldn't need to be modified
* The **config file** (e.g. `samples.yaml`)
  * Contains paths to the raw input files and other files used by the pipeline
  * We can create different configs (e.g. `biol173-multi-pop.yaml`) to specify different input files or datasets

Key Scripts:
* `update_sample_config.py` — utility script providing an easy way to create or modify config files.
* `submit_snakemake.sh` — wrapper batch file that allows us to run the pipeline in the background using SLURM.


## Running the Pipeline

### Snakemake Dry Run

Snakemake allows you to perform a **dry run** using the `-n` flag. This simply prints the steps/jobs without actually running them. It's useful for a sanity check before a real run and should be done directly on the command line as follows.

First, activate your conda environment containing Snakemake. Then run the command:
```
snakemake -np --profile {path/to/profile_directory/} --configfile {path/to/congfigfile.yaml}
```
Example usage:
```
snakemake -np --profile config/ --configfile config/samples.yaml
```

Removing the `-n` flag would run the pipeline directly on the command line, but this is inconvenient for long jobs. Instead, use SLURM.

### Snakemake Run on SLURM

To run the Snakemake pipeline in the background (as a SLURM job), use the `submit_snakemake.sh` script:
```
sbatch snakemake/submit_snakemake.sh {path/to/configfile}
```
Example usage:
```
sbatch snakemake/submit_snakemake.sh config/biol173-multi-pop.yaml
```
And that's it! To check status, use `squeue` to look for your jobs, `tail <out_file>` to check logs, or `ls <dirs_being_created>` to check outputs.

The pipeline stores outputs in a subdirectory (named per the config file) within `results/`.

View plots/results by copying them to your machine using the linux `scp` command, or view them through a virtual desktop.


## Creating a New Config
Let's say we have **new samples** in a folder `data/` that we want to process.

First we need to create a new `.yaml` to point to our samples. We can use the utility script, `update_sample_config.py`, to do this easily without editing the file directly.

The script takes a few arguments and has optional flags. To see options, run:
```
python scripts/update_sample_config.py --help
```
Example usage:
```
python scripts/update_sample_config.py data/ snakemake/samples_config/ --filename new_config
```
The script should print what it's doing to the terminal. Note that this script CANNOT change paths to the ref and anc genomes, so these must be changed manually if needed.

Once the config is successfully created, run the pipeline on it as specified above.
