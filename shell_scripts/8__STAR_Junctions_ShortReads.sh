#!/bin/bash
#SBATCH --array=1-1
#SBATCH --partition=largenode
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=80G
#SBATCH --output=star-%A_%a.out
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################

#Jenny Smith
#Run STAR aligner to provide junction files from short reads if available
#The genome index for STAR was run previously, as it is only required to create this once.

#EXAMPLE USAGE:
#ls -1 *r1.* | awk '{print $1" "$1}' | awk '{gsub("_r1","_r2",$2); print}' > samples_for_star.txt
#sbatch 8__STAR_Junctions_ShortReads.sh samples_for_star.txt

#set script to exit 1 if any of the following are not met.
set -euo pipefail

#Load modules
module purge
ml STAR/2.6.1c-foss-2016b
ml awscli/1.16.122-foss-2016b-Python-2.7.15

#Define file locations
BUCKET="s3://fh-pi-meshinchi-s"
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3/CBL/FASTQs"
cd $SCRATCH

#define Reference files
GTF="$TARGET/Reference_Data/GRCh38/gtf/gencode.v29.annotation.gtf"
GENOME="$TARGET/Reference_Data/GRCh38/fasta/genome/Gencode_GRCh38.primary_assembly.genome.fa"

#indexing commands: `STAR --runThreadN 2 --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $GENOME --sjdbGTFfile $GTF --sjdbOverhang 74`
genomeDir="$TARGET/Reference_Data/GRCh38/STAR_idx"

#Samples
sample_list=$1 #text file with 2 columns. 1 for read 1 fastq, 2 for read 2 fastq with no header.


#Define function for STAR
run_star() {

	#input fqs
	r1=$1; echo $r1
	r2=$2; echo $r2
	samp=${r1%%_*}

	#2-pass mode alignment with chimeric read detection
	#at least 25 bp of one read of a pair maps to a different loci than the rest of the read pair
	#require 20 pb on either side of chimeric junction
	#include chimeric reads in the output BAM
	#don't include chimeras that map to reference genome contig with Ns
	STAR --runMode alignReads --runThreadN 4 --genomeDir $genomeDir \
		--readFilesIn $SCRATCH/$r1 $SCRATCH/$r2 \
		--readFilesCommand zcat \
		--outFileNamePrefix ${samp}_ \
		--outSAMtype BAM SortedByCoordinate \
		--limitBAMsortRAM 63004036730 \
		--outSAMattributes NH HI NM MD AS nM jM jI XS \
		--chimSegmentMin 25 \
		--chimJunctionOverhangMin 20 \
		--chimOutType WithinBAM \
		--chimFilter banGenomicN \
		--chimOutJunctionFormat 1 \
		--twopassMode Basic \
		--twopass1readsN -1 #use all reads
}

#Select the fastq files and download from S3 Bucket
fastqs=$(cat "$sample_list" | sed -n "$SLURM_ARRAY_TASK_ID"p)
aws s3 cp --only-show-errors $BUCKET/SR/picard_fq2/$(echo "$fastqs" | awk '{print $1}' ) .
aws s3 cp --only-show-errors $BUCKET/SR/picard_fq2/$(echo "$fastqs" | awk '{print $2}' ) .

#run star aligner
run_star $fastqs

#clean up the files
rm $fastqs
