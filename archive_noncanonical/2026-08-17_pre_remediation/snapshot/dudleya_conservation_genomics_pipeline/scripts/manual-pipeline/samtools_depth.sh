#!/bin/bash
# 
#SBATCH --job-name=du_samtools_depth
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
echo "Samtools depth script started at $(date)"
#
# store params passed in as:  pipeline.sh param1 param2
#
bam=$1     # path to bam.filelist
#
# log vars
echo "Received parameter 1 (path to bam.filelist) as: $bam"
#
# Load Anaconda to run the program 
module load Anaconda3
# Activate the environment
conda activate /WAVE/projects/whittalllab/conda_envs/genomics/
echo "Conda loaded and env activated"
#
# ------- Run Samtools Depth ------- 
# FOR NOW, JUST MAKE DIR IF DNE
mkdir -pv samtools_depth
#
echo "Running depth..."
samtools depth -f $bam -q 20 > samtools_depth/depth.txt
echo "Ran depth"
#
echo "Samtools depth script finished at $(date)"