#!/bin/bash
#SBATCH --array=1-11
#SBATCH --partition=campus
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -o sqanti2-%A_%a.%j.out
#SBATCH -e sqanti2-%A_%a.%j.stderr
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################

#Jenny Smith
#Run SQANTI for annotation of isoforms from Iso-seq3 pipeline

#EXAMPLE USAGE:
# ls -1 *demux*.fastq > fastqs.txt
# ls -1 *sj.out.tab > junctions.txt
# ls -1

#set script to exit 1 if any of the following are not met.
set -euo pipefail

#Set-up environment Load modules
module purge
export PATH=~/anaconda2/bin:~/scripts/opt/bin:$PATH #for minimap2 and anaconda
export PATH=~/scripts/downloaded_software/SQANTI2:$PATH #for SQANTI2
export PATH=~/scripts/downloaded_software/cDNA_Cupcake/sequence:$PATH
export PYTHONPATH=$PATH #SQANTI2 requires PYTHONPATH to be set with cDNA_cupcake/sequence/*.py findable
ml R/3.5.3-foss-2016b-fh1
ml perl/5.22.0
source activate anaconda2.7

#Define file locations
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq"
cd $SCRATCH

#Download the CAGE peaks and splice junctions references
# curl -LO https://github.com/Magdoll/images_public/blob/master/SQANTI2_support_data/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz
curl -LO https://github.com/Magdoll/images_public/raw/master/SQANTI2_support_data/hg38.cage_peak_phase1and2combined_coord.bed.gz
curl -LO https://github.com/Magdoll/images_public/blob/master/SQANTI2_support_data/human.polyA.list.txt

#define Reference files
GTF="$TARGET/Reference_Data/GRCh38/gtf/gencode.v29.annotation.gtf"
GENOME="$TARGET/Reference_Data/GRCh38/fasta/genome/Gencode_GRCh38.primary_assembly.genome.fa"
# JUNC="$SCRATCH/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz"
CAGE="$SCRATCH/hg38.cage_peak_phase1and2combined_coord.bed.gz"
POLYA="$SCRATCH/human.polyA.list.txt"

#Define Samples
fastqs=$1 #file with all demuxed hq, collapsed, filtered fastqs files listed
juncs=$2 #STAR Aligned Junction files listed
expnfiles=$3 #Kallisto quant TPMs files listed
fl_count=$4 #e.g. test.polished.hq.no5merge.collapsed.filtered.mapped_fl_count.txt

#Function for SQANTI
sqanti_isoforms (){

	#create output directory
	samp=${hq_fq%.fastq}
	mkdir -p ${samp}

	#Run sqanti_qc.py
	sqanti_qc2.py -t 8 -z \
		--cage_peak $CAGE  \
		--polyA_motif_list $POLYA \
		--coverage $junc \
		--expression $expn \
		--fl_count $fl_count \
		--output $samp \
		--dir $SCRATCH/$samp/ $hq_fq $GTF $GENOME

	#create output directory
	mkdir -p ${samp}_filter

	#run sqanti_filter.py
	sqanti_filter.py \
	-d ${samp}_filter \
	-i $samp/${samp}_corrected.fasta $samp/${samp}_classification.txt
}

#Run the sqanti QC steps
hq_fq=$(echo "$fastqs" | sed -n "$SLURM_ARRAY_TASK_ID"p)
junc=$(echo "$juncs" | sed -n "$SLURM_ARRAY_TASK_ID"p)
expn=$(echo "$expnfiles" | sed -n "$SLURM_ARRAY_TASK_ID"p)
sqanti_isoforms $hq_fq $junc $expn $fl_count

#Clean up unncessary files
rm $CAGE
rm $POLYA
