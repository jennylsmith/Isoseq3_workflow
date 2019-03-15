#!/bin/bash
#SBATCH --partition=largenode
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=62G
#SBATCH -o lima.%j.out
#SBATCH -e lima.%j.stderr
#SBATCH -t 0-7
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################


#Jenny Smith
#Run Isoseq3 lima to remove adapters
#Combine each .ccs.bam into a consensus read set - so that lima runs on the combined .ccs.bams

#EXAMPLE USAGE:
# ccs_bams=$(ls *.ccs.bam)
# prefix="Test"
# sbatch Isoseq3_lima_primerRemoval.sh "$ccs_bams" barcodes.fasta $prefix


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
ccs_bams=$1 #can use $(ls *.ccs.bam) on command line
barcodes_fasta=$2
prefix=${3:-"ccs_combined"} #if not given, use ccs_combined as the default prefix

echo $ccs_bams
echo $prefix

#Create a combined cluster consensus sequence bams list (as an .xml file)
dataset create --type ConsensusReadSet --name $prefix ${prefix}.consensusreadset.xml $ccs_bams


#run lima on the combined .ccs bams (as listed in the .xml file). The merging happens 'under the hood'
#See pacbio "bam recipes" for more information
lima --isoseq --dump-clips --no-pbi -j 0 ${prefix}.consensusreadset.xml $barcodes_fasta ${prefix}.bam


#run dataset create on lima outputs
lima_bams=$(ls *primer_5p*3p.bam)
dataset create --type ConsensusReadSet --name ${prefix}_FL ${prefix}_FL.consensusreadset.xml $lima_bams


#Refine to remove polyA tails and concatemers
isoseq3 refine --require-polya ${prefix}_FL.consensusreadset.xml $barcodes_fasta ${prefix}.flnc.bam
