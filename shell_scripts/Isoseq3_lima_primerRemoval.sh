#!/bin/bash
#SBATCH --partition=campus
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH -o lima.%j.out
#SBATCH -e lima.%j.stderr
#SBATCH -t 0-1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################


#Jenny Smith
#Run Isoseq3 lima to remove adapters


#set-up enviornment
module purge
export PATH=~/anaconda2/bin:$PATH
source activate anaconda2.7


#Define File Locations
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
SCRATCH="/fh/scratch/delete90/meshinchi_s/jlsmith3"
dir="$SCRATCH/SMRTseq"

#Define Samples
ccs_bam=$1
primers_fasta=$2
samp=$(basename ${ccs_bam%.ccs.bam} )

echo $ccs_bam
echo $samp

#run lima
lima --isoseq --dump-clips --no-pbi -j 0 $ccs_bam $primers_fasta ${samp}_demux.bam


#run dataset create
#dataset create --type ConsensusReadSet ${samp}_combined_demux.consensusreadset.xml $(ls -1 *barcode*.bam)


#Refine
#isoseq3 refine --require-polya ${samp}_combined_demux.consensusreadset.xml $primers_fasta ${samp}_unpolished.flnc.bam
