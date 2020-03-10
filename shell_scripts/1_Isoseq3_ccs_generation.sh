#!/bin/bash
#SBATCH --array=1-25
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=62G
#SBATCH --partition=largenode
#SBATCH -o ccs.%A_%a.out
#SBATCH -e ccs.%A_%a.stderr
#SBATCH -t 0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################


#Jenny Smith
#Run Isoseq3 CCS (consensus cluster sequence)
#Each individual movie needs to have ccs.bam created.

#EXAMPLE USAGE:
# find $PWD -type f -name "*subreads.bam" > infiles.txt
# sbatch 1_Isoseq3_ccs_generation.sh infiles.txt

#set-up enviornment
module purge
source ~/scripts/Isoseq3_workflow/shell_scripts/conda.sh
conda activate isoseq3
conda env list

#set script to exit 1 if any of the following are not met.
set -euo pipefail

#Define File Locations
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3"
dir="$SCRATCH/SMRTseq/0_ccs"

#Define Samples
subread_bams=$1 #full filepaths
subread_bam=$(cat "$subread_bams" | sed -n "$SLURM_ARRAY_TASK_ID"p )
samp=$(basename ${subread_bam%.subreads.bam} ) #since given a full file path for $1, use basename

echo $subread_bam
echo $samp

#Run CCS algorithm
#Note that for isoseq3 starting version 3.2, we run Polish for CCS!
ccs --version
ccs  $subread_bam $dir/${samp}.ccs.bam --reportFile $dir/${samp}_ccs_report.txt \
	--minPasses 1 --num-threads 16 --min-rq 0.9
