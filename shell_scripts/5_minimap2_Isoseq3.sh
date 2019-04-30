#!/bin/bash
#SBATCH --partition=largenode
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=100G
#SBATCH -o minimap2.%j.out
#SBATCH -e minimap2.%j.stderr
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################

#Jenny Smith
#Minimap2 Alignment of Iso-seq3 Results for use with Cupcake Collapse Isoforms and count isoforms

#EXAMPLE USAGE:
#prefix="test"
#sbatch 5_minimap2_Isoseq3.sh  $prefix

#set script to exit 1 if any of the following are not met.
set -euo pipefail

#Load Modules/define paths. Using python2.7 in virtual enviornment "anaconda2.7"
module purge
export PATH=~/scripts/opt/bin:$PATH #for minimap2
ml SAMtools/1.8-foss-2016b
ml BamTools/2.4.1-foss-2016b

#define file locations
DELETE90="/fh/scratch/delete90/meshinchi_s/jlsmith3"
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq"
cd $SCRATCH

#define samples and genomes
genome="$TARGET/Reference_Data/GRCh38/minimap2_idx/Gencode_GRCh38.primary_assembly.genome.mmi"
prefix=${1:- "ccs_combined"} #requires same prefix that was used in 4A/B_isoseq3_polish_isoforms.sh step

#Combine the HQ FASTQs from the polish step
#this assumes you will be using all the hq.fasta/q files in the current working directory
files=$(echo ${prefix}*polished.hq.fastq*)
echo "$files"

#https://stackoverflow.com/questions/14765569/test-if-multiple-files-exist
if [[ !  $(ls -1 $files  2>/dev/null) ]]
then
	echo "No high quality fastqs in working directory."
	exit 1
elif [[ $(echo $files | wc -w) -gt 1 ]]
then
	printf "Unzipping and cating the $(echo $files | wc -w) polished hq fastqs.\n"
	cat ${prefix}*cluster_report.csv | sort | uniq > ${prefix}.polished.cluster_report.csv #sort and uniq is to remove the extra header lines "cluster_id,read_id,read_type" from the individual files
	zcat -v ${prefix}*.hq.fastq.gz > ${prefix}.polished.hq.fastq
	zcat -v ${prefix}*.lq.fastq.gz > ${prefix}.polished.lq.fastq
	ls -1 ${prefix}*.polished.bam > list
	bamtools merge -list list -out ${prefix}.polished.bam && rm list
else
	printf "unzipping the polished hq fastq.\n"
	gunzip $files
fi

#Run alignment. Fastqs can be either gzipped or not in minimap2
printf "Starting alignment.\n"
sam=${prefix}.polished.hq.fastq.sam
bam=${sam%.sam}.bam
minimap2 -ax splice -t 16 -uf --secondary=no $genome ${prefix}.polished.hq.fastq > $sam

#sort and index Bam/Sam files
printf "Sorting and Indexing.\n"
sort -k 3,3 -k 4,4n  $sam > ${sam%.sam}.srt.sam
samtools view -bS $sam | samtools sort -o ${bam%.bam}.srt.bam - && rm $sam
samtools index ${bam%.bam}.srt.bam
