#!/bin/bash
#SBATCH --partition=largenode
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=62G
#SBATCH -o minimap2.%j.out
#SBATCH -e minimap2.%j.stderr
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################

#Jenny Smith
#Minimap2 Alignment of Iso-seq3 Results for use with Cupcake Collapse Isoforms and count isoforms


#Load Modules/define paths. Using python2.7 in virtual enviornment "anaconda2.7"
module purge
source activate anaconda2.7
export PATH=~/scripts/opt/bin:$PATH #for minimap2 & cDNA_cupcake
ml SAMtools/1.8-foss-2016b


#define file locations
DELETE90="/fh/scratch/delete90/meshinchi_s/jlsmith3"
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq"
cd $scratch

#define samples and genomes
genome="$TARGET/Reference_Data/GRCh38/fasta/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
prefix=${1:- "minimap2"}


#Combine the HQ FASTQ from the polish step
#this assumes you will be using all the hq.fasta/q files in the current working directory
# and requires a prefix that was used in 4_isoseq3_polish_isoforms.sh step.
echo ${prefix}*polished.hq.fastq.gz
files=$(echo ${prefix}*polished.hq.fastq.gz)

# if [[ ! -e "$files" ]]
# then
# 	echo "No high quality fastqs in working directory."
# 	exit 1
# else if [[ $(echo $files | wc -l) -gt 1 ]]
	gunzip ${prefix}*polished.hq.fastq.gz)
	cat ${prefix}*.hq.fastq > ${prefix}.polished.hq.fastq
# else
# 	gunzip ${prefix}.polished.hq.fastq.gz
# fi


#Run alignment
minimap2 -ax splice -t 16 -uf --secondary=no $genome ${prefix}_polished.hq.fastq > ${prefix}_polished.hq.fastq.sam

#sort and index Bam/Sam files
sort -k 3,3 -k 4,4n  ${prefix}_polished.hq.fastq.sam > ${prefix}_polished.hq.fastq.srt.sam
# samtools view -bS ${hq_fq}.sam > ${hq_fq}.bam
# samtools sort ${hq_fq}.bam > ${hq_fq}.srt.bam
# samtools index ${hq_fq}.srt.bam
