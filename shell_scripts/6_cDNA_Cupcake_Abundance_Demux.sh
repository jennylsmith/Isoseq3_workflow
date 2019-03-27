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
source activate anaconda2.7
export PATH=~/.local/bin:$PATH #for cDNA_cupcake
ml SAMtools/1.8-foss-2016b


#define file locations
DELETE90="/fh/scratch/delete90/meshinchi_s/jlsmith3"
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq"
manifest=
cd $scratch

#define samples
sample_idmap="~/scripts/Isoseq3_workflow/samples/Sample_ID_Map.csv" #location from manifest.sh
hq_fq=$1
prefix=${2:-"ccs_combined"} #prefix must be the same as used in 44A_isoseq3_polish_isoforms.sh


#Collapse Reads
collapse_isoforms_by_sam.py --input ${hq_fq} --fq \
   -s ${hq_fq}.srt.sam -c 0.99 -i 0.95 --dun-merge-5-shorter \
   -o ${prefix}.polished.hq.no5merge

#calculate abundance
get_abundance_post_collapse.py ${prefix}.polished.hq.no5merge.collapsed  ${prefix}.polished.cluster_report.csv

#Filter away 5' degraded isoforms
filter_away_subset.py ${prefix}.polished.hq.no5merge.collapsed





#Run demultiplexing scripts
python <path_to_cupcake>/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fastq test_polished.hq.no5merge.collapsed.filtered.rep.fq --read_stat test_polished.hq.no5merge.collapsed.read_stat.txt --classify_csv ../flnc.report.hacked.csv -o test_polished.hq.no5merge.collapsed.filtered.mapped_fl_count.txt
