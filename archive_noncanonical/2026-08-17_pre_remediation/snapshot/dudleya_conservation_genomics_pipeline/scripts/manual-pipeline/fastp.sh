#!/bin/bash
# 
#SBATCH --job-name=du_fastp
#SBATCH --output=out-fastp-%j.out 
# 
#SBATCH --partition=cmp 
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16 
#SBATCH --mem-per-cpu=1G 
#SBATCH --time=2:00:00 
#
#SBATCH --mail-user=ehackstadt@scu.edu
#SBATCH --mail-type=END
#
set -e
echo "Fastp script started at $(date)"
#
# store params passed in as:  pipeline.sh param1 param2
#
read1=$1     # path to R1 fastq.gz file
read2=$2     # path to R2 fastq.gz file
dir=$3      # path to working directory to save outputs
#
# utility variables for filenames
r1="_R1"
r2="_R2"
name=$(basename $read1 _R1.fastq.gz)
# log vars
echo "Received parameter 1 (path to read 1) as: $read1"
echo "Recieved parameter 2 (path to read 2) as: $read2"
echo "Recieved parameter 3 (working directory) as: $dir"
echo "Parsed sample name as: $name"
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
# ------- Run Fastp ------- 
#
# clear old files — DISABLED FOR NOW
# if [ -d "fastp_filtering/" ]; then
#     rm -rv fastp_filtering/
# fi
# mkdir -v fastp_filtering
#
# FOR NOW, JUST MAKE DIR IF DNE
mkdir -pv fastp_filtering
#
echo "Running fastp..."
fastp \
-i $read1 \
-I $read2 \
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
echo "Fastp script finished at $(date)"