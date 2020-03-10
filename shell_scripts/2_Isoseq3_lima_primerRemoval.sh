#!/bin/bash
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=62G
#SBATCH --partition=largenode
#SBATCH --mem=62G
#SBATCH -o lima.%j.out
#SBATCH -e lima.%j.stderr
#SBATCH -t 0-7
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################


#Jenny Smith
#Run Isoseq3 lima to remove adapters
#Combine each .ccs.bam into a consensus read set - so that lima runs on the combined .ccs.bams

#EXAMPLE USAGE:
# ccs_bams=$(find $PWD -name "*.ccs.bam") #full file paths required
# prefix="Test"
# sbatch 2_Isoseq3_lima_primerRemoval.sh "$ccs_bams" barcodes.fasta $prefix


#set-up enviornment
module purge
source ~/scripts/Isoseq3_workflow/shell_scripts/conda.sh
conda activate isoseq3
conda env list
ml BamTools/2.4.1-foss-2016b


#set script to exit 1 if any of the following are not met.
#Still cannot set at beginning of script due to conda error $PS1 variable unbound...
set -euo pipefail

#Define File Locations
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3"
dir="$SCRATCH/SMRTseq/1_lima"

#Define Samples
ccs_bams=$1 #can use $(ls *.ccs.bam) on command line
barcodes_fasta=$2
prefix=${3:-"ccs_combined"} #if not given, use ccs_combined as the default prefix

echo $ccs_bams
echo $prefix
echo $barcodes_fasta

#Create a combined cluster consensus sequence bams list (as an .xml file)
if [ ! -e $dir/${prefix}.ccs.consensusreadset.xml ]
then
dataset create --type ConsensusReadSet --name $prefix $dir/${prefix}.ccs.consensusreadset.xml $ccs_bams
fi

#run lima on the combined .ccs bams (as listed in the .xml file). The merging happens 'under the hood'
#Use --peek-guess to remove spurious matches (only applicable if you supply multiple primer pairs).
#See pacbio "bam recipes" for more information
lima --version
if ! ls *primer_5p*3p.bam 2>&1 /dev/null
then
lima --isoseq --dump-clips --peek-guess -j 16 --log-file ${prefix}.lima.log \
	$dir/${prefix}.ccs.consensusreadset.xml $barcodes_fasta $dir/${prefix}.bam
fi
echo "lima complete"


#run dataset create on lima outputs.
#Dataset now runs only in python3 environement 
cd $dir
if [ ! -e ${prefix}_FL.consensusreadset.xml ]
then
conda activate pacbio_py3
lima_bams=$(ls *primer_5p*3p.bam)
dataset create --force --type ConsensusReadSet --name ${prefix}_FL ${prefix}_FL.consensusreadset.xml $lima_bams
conda deactivate
fi
echo "Done creating consensus read set"

#Refine to remove polyA tails and concatemers
isoseq3 refine -j 16 --log-file ${prefix}.flnc.log --require-polya \
	${prefix}_FL.consensusreadset.xml $barcodes_fasta ${prefix}.flnc.bam
echo "Isoseq3 refine complete"

#Convert to fastq
bamtools convert -format fastq -in ${prefix}.flnc.bam > ${prefix}.flnc.fastq
echo "conversion of FLNC.bam to fastq complete"
