#!/bin/bash
#SBATCH --partition=largenode
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=62G
#SBATCH -o clust.%j.out
#SBATCH -e clust.%j.stderr
#SBATCH -t 0-1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################


#Jenny Smith
#Run Isoseq3 Cluster to cluster full-length non-concatemer reads into isoforms

#EXAMPLE USAGE:


#set script to exit 1 if any of the following are not met.
set -euo pipefail


#set-up enviornment
module purge
export PATH=~/anaconda2/bin:$PATH
source activate anaconda2.7


#Define File Locations
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3"
dir="$SCRATCH/SMRTseq"


#Define Samples
flnc_bam=$1
prefix=${2:-"ccs_combined"}

echo $flnc_bam
echo $prefix

#run isoseq3 cluster step
isoseq3 cluster -j 16 --verbose --log-file "${prefix}_cluster.log" --split-bam 24 $flnc_bam ${prefix}.unpolished.bam
