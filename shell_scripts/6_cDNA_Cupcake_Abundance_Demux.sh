#!/bin/bash
#SBATCH --partition=largenode
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=62G
#SBATCH -o cupcake.%j.out
#SBATCH -e cupcake.%j.stderr
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################

#Jenny Smith
#cDNA Cupcake Collapse Isoforms and count isoforms from Isoseq3 workflow


#Load Modules/define paths. Using python2.7 in virtual enviornment "anaconda2.7"
module purge
export PATH=~/anaconda2/bin:~/scripts/Isoseq3_workflow/shell_scripts:~/scripts/downloaded_software/cDNA_Cupcake/post_isoseq_cluster:$PATH #for ananconda, cDNA_cupcake, and custom flnc_report.py
source activate anaconda2.7

#define file locations
# CDNA_CUPCAKE="~/scripts/downloaded_software/cDNA_Cupcake"
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq/testing_set"
cd $SCRATCH

#define reference files
manifest="/home/jlsmith3/scripts/Isoseq3_workflow/samples/Sample_ID_Map.csv" #same location from 00_manifest.sh

#define samples
hq_fq=$1 #eg hq_fq=test.polished.hq.fastq
prefix=${2:-"ccs_combined"} #same prefix as used in  3_Isoseq3_cluster_isoforms.sh, 4A/B_isoseq3_polish_isoforms.sh , and 5_minimap2_Isoseq3.sh

#Collapse Reads
collapse_isoforms_by_sam.py --input ${hq_fq} --fq \
   -s ${hq_fq}.srt.sam -c 0.99 -i 0.95 --dun-merge-5-shorter \
   -o ${prefix}.polished.hq.no5merge

#calculate abundance
get_abundance_post_collapse.py ${prefix}.polished.hq.no5merge.collapsed  ${prefix}.polished.cluster_report.csv

#Filter away 5' degraded isoforms
filter_away_subset.py ${prefix}.polished.hq.no5merge.collapsed

#create collapsed flnc.report.csv where each pair of movies is defined its sample condition, and no longer by movie ID (eg m54228_181211_220100)
#This assumes flnc.report.csv is in the current working directoy.
6__hack_flnc_report.py "$PWD/${prefix}.flnc.report.csv" $manifest $prefix

#Run demultiplexing scripts
classify_csv="$PWD/${prefix}.flnc.report.hacked.csv"
demux_isoseq_with_genome.py \
  --mapped_fastq ${prefix}.polished.hq.no5merge.collapsed.filtered.rep.fq \
  --read_stat ${prefix}.polished.hq.no5merge.collapsed.read_stat.txt \
  --classify_csv $classify_csv \
  -o ${prefix}.polished.hq.no5merge.collapsed.filtered.mapped_fl_count.txt

#Create demultiplexed GFF and FASTA/FASTQ files
tuple=cat $manifest | cut -f 4 -d "," | grep -v "Sample_Number" | awk -v q="'" -v p="(" -v p2=")" '{print p q $1q","q$1 q p2}' | tr "\n"  "," | sed -E 's/^(.+),$/"\1"\n/'
# demux_by_barcode_groups.py --pooled_fastx ${hq_fq} ${hq_fq}.srt.sam ${prefix}.polished.hq.no5merge.collapsed.filtered.mapped_fl_count.txt ${prefix}_demux  $tuple #"('AML','AML'),('NBM','NBM')"
