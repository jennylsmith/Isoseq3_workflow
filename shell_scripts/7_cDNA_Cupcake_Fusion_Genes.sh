#!/bin/bash
#SBATCH --partition=campus
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH -o fusion.%j.out
#SBATCH -e fusion.%j.stderr
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################

#Jenny Smith
#cDNA Cupcake Collapse Isoforms and count isoforms from Isoseq3 workflow

#EXAMPLE USAGE:


#set script to exit 1 if any of the following are not met.
set -euo pipefail

#Load Modules/define paths. Using python2.7 in virtual enviornment "anaconda2.7"
module purge
export PATH=~/anaconda2/bin:~/scripts/Isoseq3_workflow/shell_scripts:$PATH #for ananconda, and custom flnc_report.py
export PATH=~/scripts/downloaded_software/cDNA_Cupcake/post_isoseq_cluster:~/scripts/opt/bin:$PATH
ml SAMtools/1.9-foss-2016b
ml picard/2.18.29-Java
source activate anaconda2.7

#define file locations
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq"
cd $SCRATCH

#define samples
hq_fq=$1 #eg test.polished.hq.fastq
prefix=${2:-"ccs_combined"} #same prefix as used in  3_Isoseq3_cluster_isoforms.sh, 4A/B_isoseq3_polish_isoforms.sh , and 5_minimap2_Isoseq3.sh

# #Run fusion detection algorithm
# printf "Running fusion dection step.\n"
# fusion_finder.py --input $hq_fq --fq -s ${hq_fq}.srt.sam \
#   -o ${prefix}_isoforms.fasta.fusion \
#   --cluster_report_csv ${prefix}.polished.cluster_report.csv

#Filter away 5' degraded isoforms
# filter_away_subset.py ${prefix}_isoforms.fasta.fusion

#Run demultiplexing scripts
# --classify_csv file created from 6_cDNA_Cupcake_Abundance_Demux.sh
printf "Demultiplexing by sample condition.\n"
demux_isoseq_no_genome.py --hq_fastq ${prefix}_isoforms.fasta.fusion.rep.fq \
    --cluster_csv ${prefix}.polished.cluster_report.csv \
    --classify_csv ${prefix}.flnc.report.hacked.csv \
    -o ${prefix}_isoforms.fusion.fl_count.txt

#Use minimap2 to create a pooled sam file from the pooled fastq
sam=${prefix}_isoforms.fasta.fusion.rep.fq.sam
java -jar $EBROOTPICARD/picard.jar FastqToSam FASTQ=${prefix}_isoforms.fasta.fusion.rep.fq O=$sam READ_GROUP_NAME=PBfusion SAMPLE_NAME=Pooled PLATFORM=PacBio

#Create demultiplexed GFF and FASTA/FASTQ files
printf "Creating sample level GFF, Fastqs, and BAM files.\n"
tuple=$(cat $manifest | cut -f 3 -d "," | grep -v "Reg." | sort | uniq | awk -v q="'" -v p="(" -v p2=")" '{print p q $1 q "," q $1 q p2}' | tr "\n"  "," | sed -E 's/^(.+),$/"\1"\n/')
demux_by_barcode_groups_fusion.py --pooled_fastx ${sam%.sam} $sam ${prefix}_isoforms.fusion.fl_count.txt ${prefix}_fusion_demux  $tuple

#align, Sort and index the demuxed individual files
# for file in $(ls -1 ${prefix}_demux*fastq)
# do
#      printf "Starting to map $file."
#      minimap2 -ax splice -t 4 -uf --secondary=no $genome $file | samtools view -bS - | samtools sort -o ${file%.fastq}.srt.bam -
#      samtools index ${file%.fastq}.srt.bam
# done
