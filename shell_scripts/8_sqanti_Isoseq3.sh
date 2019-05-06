#!/bin/bash
#SBATCH --array=1-11
#SBATCH --partition=campus
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH -o sqanti2-%A_%a.out
#SBATCH -e sqanti2-%A_%a.stderr
#SBATCH -t 7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################

#Jenny Smith
#Run SQANTI for annotation of isoforms from Iso-seq3 pipeline

#EXAMPLE USAGE:
# paste -d " " <(ls -1 *.fastq)  <(ls -1 *SJ.out.tab) > sqanti_files.txt
# sbatch 8_sqanti_isoseq3.sh sqanti_files.tx test.polished.hq.no5merge.collapsed.filtered.mapped_fl_count.txt

#set script to exit 1 if any of the following are not met.
set -euo pipefail

#Set-up environment Load modules
module purge
export PATH=~/anaconda2/bin:~/scripts/opt/bin:$PATH #for minimap2 and anaconda
export PATH=~/scripts/downloaded_software/SQANTI2:$PATH #for SQANTI2
export PATH=~/scripts/downloaded_software/cDNA_Cupcake/sequence:$PATH #for cDNA_cupcake
export PYTHONPATH=$PATH #SQANTI2 requires PYTHONPATH to be set with cDNA_cupcake/sequence/*.py findable
ml R/3.5.3-foss-2016b-fh1 #for final QC plots
ml perl/5.22.0 #for geneMarkS program
source activate anaconda2.7

#Define file locations
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq"
cd $SCRATCH

#Need only Once: Download the CAGE peaks, ployA, and intropolis junctions references
# curl -L https://github.com/Magdoll/images_public/raw/master/SQANTI2_support_data/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.gz | gunzip -c > intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.bed
# curl -L https://github.com/Magdoll/images_public/raw/master/SQANTI2_support_data/hg38.cage_peak_phase1and2combined_coord.bed.gz | gunzip -c > hg38.cage_peak_phase1and2combined_coord.bed
# curl -LO https://raw.githubusercontent.com/Magdoll/images_public/master/SQANTI2_support_data/human.polyA.list.txt

#Define Reference files
GTF="$TARGET/Reference_Data/GRCh38/gtf/gencode.v29.annotation.gtf"
GENOME="$TARGET/Reference_Data/GRCh38/fasta/genome/Gencode_GRCh38.primary_assembly.genome.fa"
CAGE="$SCRATCH/hg38.cage_peak_phase1and2combined_coord.bed"
POLYA="$SCRATCH/human.polyA.list.txt"
JUNCS="$SCRATCH/intropolis.v1.hg19_with_liftover_to_hg38.tsv.min_count_10.modified.bed"

#Define Samples
sample_file_names=$1 #a 2 column file, space delimited. no header. Col1 = demux.fastq, Col2 = SJ.out.tab)
fl_count=$2 #e.g. test.polished.hq.no5merge.collapsed.filtered.mapped_fl_count.txt

#Function for SQANTI
sqanti_isoforms (){

	#define files
	fq=$1; junc=$2; fl=$3
	echo $fq

	#create output directory and junctions directory
	samp=${fq%.fastq}
	mkdir -p $SCRATCH/${samp}/junction_refs
	cp -n $junc $SCRATCH/${samp}/junction_refs/${junc}.bed #rename to SJ.tab to bed file
	ln -s $JUNCS $SCRATCH/${samp}/junction_refs #symlink to intropolis reference

	#Run sqanti_qc.py
	sqanti_qc2.py $fq $GTF $GENOME -t 8 \
		--cage_peak $CAGE  \
		--polyA_motif_list $POLYA \
		--coverage "$SCRATCH/${samp}/junction_refs/*.bed" \
		--fl_count $fl \
		--output ${samp} \
		--dir $SCRATCH/$samp/ 2> $SCRATCH/${samp}/${samp}_qc.log
		# --expression $3 \ #not implemented yet in sqanti2

	#Ensure filter results go in the same output directory
	cd $SCRATCH/${samp}/

	#run sqanti_filter.py (missing interpreter in pyhton script, hence full path required)
	python $(which sqanti_filter2.py) -a 0.8 -c 1 \
		$SCRATCH/$samp/${samp}_classification.txt \
		$SCRATCH/$samp/${samp}.renamed_corrected.fasta \
		$SCRATCH/$samp/${samp}.renamed_corrected.sam 2> $SCRATCH/${samp}/${samp}_filter.log

}

#Must divide the demultiplexed count file by sample
infiles=$(cat $sample_file_names | sed -n "$SLURM_ARRAY_TASK_ID"p)
demux_sample_cond=$(echo $infiles | sed -E 's/^.+demux_(.+)_only.+/\1/')
col_number=$(head -1 $fl_count | tr "," "\n" | sed -n "/$demux_sample_cond/=")
filename=${infiles%%.fastq*}_collapsed.filtered.mapped_fl_count.txt
cat $fl_count | sed -E 's/id/pbid/' | cut -f 1,$col_number -d "," | tr "," "\t" > $filename

#Run the sqanti QC steps
sqanti_isoforms $infiles $filename
