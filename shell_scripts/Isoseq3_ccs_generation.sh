#!/bin/bash
#SBATCH --partition=largenode
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=50G
#SBATCH -o ccs.%j.out
#SBATCH -e ccs.%j.stderr
#SBATCH -t 0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################



#Jenny Smith
#Run Isoseq3 CCS (consensus cluster sequence)

#set-up enviornment
module purge
export PATH=~/anaconda2/bin:$PATH
source activate anaconda2.7



#Define File Locations
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3"
dir="$SCRATCH/SMRTseq"

#Define Samples
subread_bam=$1
samp=$(basename ${subread_bam%.subreads.bam} ) #if given a full file path for $1 use basename

echo $subread_bam
echo $samp
#Run CCS algorithm
# --reportFile "${samp}.ccs_report.txt" #reportFile is not recognized? yet all mine got overwritten....

ccs  $subread_bam ${samp}.ccs.bam --reportFile ${samp}_ccs_report.txt --noPolish --minPasses 1 --numThreads 16
