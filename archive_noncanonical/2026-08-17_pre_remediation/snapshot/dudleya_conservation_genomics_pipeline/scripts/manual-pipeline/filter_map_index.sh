#!/bin/bash
# 
#SBATCH --job-name=du_pipeline
#SBATCH --output=out-pipeline-%j.out 
# 
#SBATCH --partition=condo 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16 
#SBATCH --mem-per-cpu=1G
#SBATCH --time=4:00:00 
#
#SBATCH --mail-user=ehackstadt@scu.edu
#SBATCH --mail-type=END
#
set -e
echo "Pipeline script started at $(date)"
#
# store params passed in as:  pipeline.sh param1 param2
#
name=$1     # string representing a custom name, e.g. DU059LP015 or DU229
dir=$2      # path to working directory containing R1 and R2 fastq.gz files
# path to the reference genome .fasta file. if non provided, uses:
ref=${3:-'/WAVE/projects/whittalllab/dudleya/draft_genome/Dudleya_hifiasm_purged_host_genome.fasta'}
#
# utility variables for filenames
r1="_R1"
r2="_R2"
# log vars
echo "Received parameter 1 (name for outputs) as: $name"
echo "Recieved parameter 2 (working dir) as: $dir"
echo "Recieved parameter 3 (path to ref genome) as: $ref"
#
# Load Anaconda to run the program 
module load Anaconda3
# Activate the environment
conda activate /WAVE/projects/whittalllab/conda_envs/genomics/
echo "Conda loaded and env activated"
# Folder structure & clear old files
cd $dir
echo "Working dir initially contains:"
ls
#
# ------- Run the pipeline commands ------- 
#
# Fastp
# clear old files
if [ -d "fastp_filtering/" ]; then
    rm -rv fastp_filtering/
fi
mkdir -v fastp_filtering
echo "Running fastp..."
fastp \
-i *R1*.fastq*.gz \
-I *R2*.fastq*.gz \
-o fastp_filtering/filtered-polyG_$name$r1.fastq.gz \
-O fastp_filtering/filtered-polyG_$name$r2.fastq.gz \
--trim_poly_g \
--length_required 100 \
-q 20 \
--thread 4 \
--html fastp_filtering/polyG_FastpFilterReport_$name.html \
--json fastp_filtering/polyG_FastpFilterReport_$name.json
ls -lh fastp_filtering
echo "Ran fastp"
#
# BWA
# clear old files
if [ -d "aln/" ]; then
    rm -rv aln/
fi
mkdir -v aln
echo "Running BWA..."
# run BWA
# first check if ref genome is already indexed; if not, index it
# if [ ! -f *.amb ]; then
#     bwa index $ref
# fi
bwa mem \
-t 8 \
$ref \
fastp_filtering/filtered-polyG_$name$r1.fastq.gz \
fastp_filtering/filtered-polyG_$name$r2.fastq.gz \
 | samtools sort -o aln/$name-sorted.bam
samtools index aln/$name-sorted.bam
samtools flagstat aln/$name-sorted.bam > aln/$name-flagstat.txt
ls -lh aln
# manually save the commands used since BWA does not
echo "bwa mem \
-t 8 \
$ref \
fastp_filtering/filtered-polyG_$name$r1.fastq.gz \
fastp_filtering/filtered-polyG_$name$r2.fastq.gz \
 | samtools sort -o aln/$name-sorted.bam
samtools index aln/$name-sorted.bam
samtools flagstat aln/$name-sorted.bam > aln/$name-flagstat.txt" > aln/command.txt
echo "Ran BWA"
echo "Pipeline script finished at $(date)"