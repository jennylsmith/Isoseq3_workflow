#!/bin/bash
#SBATCH --array=1-24
#SBATCH --partition=largenode
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=250G
#SBATCH --output=polish-%A_%a.out
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################

#Jenny Smith
#Run Isoseq3 polish

#EXAMPLE USAGE:
# sbatch 4A_isoseq3_polish_isoforms.sh

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

#Change working directory to where the subread.bams and the unpolished.bam are located.
cd $dir
files=$(ls -1 *unpolished*bam) #this assumes all unpolished.bam files in the current working directory are included.
echo "$files"

#function to polish each file.
polish (){

	#Define Samples
	unpolished_bam=$1
	prefix=$(echo $unpolished_bam | sed -E "s/\.bam|\.unpolished//gm" ) #Uses same prefix as in 3_Isoseq3_cluster_isoforms.py
	echo $unpolished_bam $prefix

	#Create a subread set for the original samples' "movie.subreads.bam" files. Uses `find` command to create a list of the input subreads.bam files with full file paths.
	#this assumes all subreads.bam files in the current working directory are included.
	subread_bams=$(find $PWD -type f -name "*.subreads.bam"); echo "$subread_bams"
	dataset create --force --type SubreadSet --name $prefix  ${prefix}.subreadset.xml $subread_bams

	# run the serial polishing algorithm
	isoseq3 polish -j 16  --verbose --log-file "${prefix}_polish.log"  $unpolished_bam ${prefix}.subreadset.xml ${prefix}.polished.bam
}

#Run the polish steps
unpolished_bam=$(echo "$files" | sed -n "$SLURM_ARRAY_TASK_ID"p)
polish $unpolished_bam
