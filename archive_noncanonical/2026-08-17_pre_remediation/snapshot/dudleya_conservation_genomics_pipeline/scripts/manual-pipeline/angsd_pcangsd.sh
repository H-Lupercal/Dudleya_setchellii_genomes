#!/bin/bash
# 
#SBATCH --job-name=du_angsd
#SBATCH --output=out-angsd-%j.out 
# 
#SBATCH --partition=mem
#SBATCH --nodes=1 
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=16 
#SBATCH --mem-per-cpu=6G
#SBATCH --time=4:00:00 
#
#SBATCH --mail-user=ehackstadt@scu.edu
#SBATCH --mail-type=END
#
set -e
echo "ANGSD script started at $(date)"
#
# store params passed in as:  pipeline.sh param1 param2 param3 param4
#
b=$1     # path to bam.filelist
dir=$2      # path to working directory from which the bam.filelist is relative to
# path to the reference genome .fasta file. if non provided, uses:
ref=${3:-'/WAVE/projects/whittalllab/dudleya/draft_genome/Dudleya_hifiasm_purged_host_genome.fasta'}
# path to .fasta file of the ancestral genome. if none provided, uses:
anc=${4:-"/WAVE/projects/whittalllab/dudleya/ancestral_genome_Du-14/D_cymosa_consensus.fa.gz"}
#
# log vars
echo "Received parameter 1 (path to bam.filelist) as: $b"
echo "Recieved parameter 2 (path to output dir) as: $dir"
echo "Recieved parameter 3 (path to ref genome) as: $ref"
echo "Recieved parameter 4 (path to anc genome) as: $anc"
#
# Load Anaconda to run the program 
module load Anaconda3
# Activate the environment
conda activate /WAVE/projects/whittalllab/conda_envs/genomics/
echo "Conda loaded and env activated"
#
# ------- Run ANGSD & PCAngsd ------- 
#
# ANGSD
# clear old files
if [ -d "$dir/angsd_out" ]; then
    rm -rv $dir/angsd_out/
fi
mkdir -v $dir/angsd_out/
#
angsd \
-b $b \
-GL 1 \
-ref $ref \
-anc $anc \
-doGlf 2 \
-doMajorMinor 1 \
-doCounts 1 \
-doMaf 1 \
-doSaf 1 \
-minMapQ 30 \
-minQ 20 \
-doPost 1 \
-doGeno 32 \
-uniqueOnly 1 \
-SNP_pval 1e-6 \
-trim 5 \
-nThreads 8 \
-out $dir/angsd_out/angsd_out
echo "Ran ANGSD"
#
# PCAngsd
#
# clear old files
if [ -d "$dir/pcangsd_out/" ]; then
    rm -rv pcangsd_out/
fi
mkdir -v $dir/pcangsd_out/
#
pcangsd \
-b $dir/angsd_out/angsd_out.beagle.gz \
-o $dir/pcangsd_out/pcangsd_out \
--selection
echo "Ran PCAngsd"
echo "ANGSD script finished at $(date)"

