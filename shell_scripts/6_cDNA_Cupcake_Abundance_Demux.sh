#!/bin/bash
#SBATCH --partition=campus
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH -o cupcake.%j.out
#SBATCH -e cupcake.%j.stderr
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################

#Jenny Smith
#cDNA Cupcake Collapse Isoforms and count isoforms from Isoseq3 workflow and demultiplex the pooled Samples.


#EXAMPLE USAGE
#hq_fq="$PWD/test.polished.hq.fastq"
#flnc_report="$PWD/test.flnc.report.csv"
#prefix="test"
#sbatch 6_cDNA_Cupcake_Abundance_Demux.sh $hq_fq $flnc_report $prefix

#set script to exit 1 if any of the following are not met.
set -euo pipefail

#Load Modules/define paths. Using python2.7 in virtual enviornment "anaconda2.7"
module purge
export PATH=~/anaconda2/bin:~/scripts/Isoseq3_workflow/shell_scripts:$PATH #for ananconda, and custom flnc_report.py
export PATH=~/scripts/downloaded_software/cDNA_Cupcake/post_isoseq_cluster:~/scripts/opt/bin:$PATH # for cDNA_cupcake and minimap2
ml SAMtools/1.9-foss-2016b
source activate anaconda2.7

#define file locations
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq"
cd $SCRATCH

#define reference files
genome="$TARGET/Reference_Data/GRCh38/minimap2_idx/Gencode_GRCh38.primary_assembly.genome.mmi"
manifest="/home/jlsmith3/scripts/Isoseq3_workflow/samples/Sample_ID_Map.csv" #same location from 00_manifest.sh

#define samples (NTE:requires full file paths)
hq_fq=$1 #eg $PWD/test.polished.hq.fastq
flnc_report=$2 #eg $PWD/test.flnc.report.csv
prefix=${3:-"ccs_combined"} #same prefix as used in  3_Isoseq3_cluster_isoforms.sh, 4A/B_isoseq3_polish_isoforms.sh , and 5_minimap2_Isoseq3.sh

#Collapse Reads
printf "Starting to collapse isoforms.\n"
collapse_isoforms_by_sam.py --input ${hq_fq} --fq -s ${hq_fq}.srt.sam -c 0.99 -i 0.95 --dun-merge-5-shorter -o ${prefix}.polished.hq.no5merge

#calculate abundance
printf "Calculating isoform abundance.\n"
get_abundance_post_collapse.py ${prefix}.polished.hq.no5merge.collapsed ${prefix}.polished.cluster_report.csv

#Filter away 5' degraded isoforms
filter_away_subset.py ${prefix}.polished.hq.no5merge.collapsed

#create collapsed flnc.report.csv where each pair of movies is defined its sample condition, and no longer by movie ID (eg m54228_181211_220100)
6__hack_flnc_report.py $flnc_report $manifest $prefix

#Run demultiplexing scripts
printf "Demultiplexing by sample condition.\n"
demux_isoseq_with_genome.py \
  --mapped_fastq ${prefix}.polished.hq.no5merge.collapsed.filtered.rep.fq \
  --read_stat ${prefix}.polished.hq.no5merge.collapsed.read_stat.txt \
  --classify_csv ${prefix}.flnc.report.hacked.csv \
  -o ${prefix}.polished.hq.no5merge.collapsed.filtered.mapped_fl_count.txt

# # #Use minimap2 to create a pooled sam file from the pooled fastq
sam=${prefix}.polished.hq.no5merge.collapsed.filtered.rep.fq.sam
minimap2 -ax splice -t 4 -uf --secondary=no $genome ${prefix}.polished.hq.no5merge.collapsed.filtered.rep.fq > $sam

#Create demultiplexed GFF and FASTA/FASTQ files
printf "Creating sample level GFF, Fastqs, and BAM files.\n"
tuple=$(cat $manifest | cut -f 3 -d "," | grep -v "Reg." | sort | uniq | awk -v q="'" -v p="(" -v p2=")" '{print p q $1 q "," q $1 q p2}' | tr "\n"  "," | sed -E 's/^(.+),$/"\1"\n/')
demux_by_barcode_groups.py --pooled_fastx ${sam%.sam} $sam ${prefix}.polished.hq.no5merge.collapsed.filtered.mapped_fl_count.txt ${prefix}_demux  $tuple

#align, Sort and index the demuxed individual files
for file in $(ls -1 ${prefix}_demux*fastq)
do
   printf "Starting to align $file.\n"
   minimap2 -ax splice -t 4 -uf --secondary=no $genome $file | samtools view -bS - | samtools sort -o ${file%.fastq}.srt.bam -
   samtools index ${file%.fastq}.srt.bam
done

# sort -k 3,3 -k 4,4n  ${file%.fastq}.sam > ${file%.fastq}.srt.sam #&& rm ${file%.fastq}.sam
# samtools view -bS ${file%.fastq}.sam | samtools sort -o ${file%.fastq}.srt.bam -
