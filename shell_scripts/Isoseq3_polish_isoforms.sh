#!/bin/bash
#SBATCH --partition=largenode
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=62G
#SBATCH -o polish.%j.out
#SBATCH -e polish.%j.stderr
#SBATCH -t 0-1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################


#Jenny Smith
#Run Isoseq3 polish

#EXAMPLE USAGE:


#set-up enviornment
module purge
export PATH=~/anaconda2/bin:$PATH
source activate anaconda2.7


#Define File Locations
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3"
dir="$SCRATCH/SMRTseq"


#Change working directory to where the subread.bams and the unpolished.bam are located.
#subread.bams can be nested in directories - will use recursive 'find' command to get full file paths.
cd $dir


#Define Samples
unpolished_bam=$1
prefix=$2

#If prefix argument is blank on the command line - pick a default prefix name
if [ $(echo $prefix | wc -w) -eq 0 ]
then
        prefix="ccs_combined"
fi

echo $unpolished_bam

#Create a subread set for the original samples' "movie.subreads.bam" files.
#Use find command to create a list of the input subreads.bam files with full file paths.
#this assumes all subreads.bam files in the current working directory are included.
subread_bams=$(find $PWD -type f -name "*.subreads.bam")
dataset create --type SubreadSet --name $prefix ${prefix}.subreadset.xml $subread_bams

#run the serial polishing algorithm
isoseq3 polish -j 16 --verbose --log-file "${prefix}_polish.log"  $unpolished_bam ${prefix}.subreadset.xml ${prefix}_polished.bam
