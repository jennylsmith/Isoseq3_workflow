#!/bin/bash
#SBATCH --partition=campus
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -o sqanti.%j.out
#SBATCH -e sqanti.%j.stderr
#SBATCH -t 0-1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################



#Jenny Smith
#Run SQANTI for annotation of isoforms from Iso-seq3 pipeline


#Load modules
module purge
export PATH=~/scripts/opt/bin:$PATH #for GMAP & sqanti
ml Perl/5.28.0-foss-2016b
ml R/3.5.1-foss-2016b-fh1
module unload cURL/7.49.1-foss-2016b #automatically loaded with one of the other modules - but causes "ImportError: /app/easybuild/software/cURL/7.49.1-foss-2016b/lib/libcurl.so.4: version `CURL_OPENSSL_3' not found (required by /home/jlsmith3/.local/lib/python2.7/site-packages/pysam/libchtslib.so)"


#define working directory.
TARGET="/fh/fast/meshinchi_s/workingDir/TARGET"
GMAPDB="/fh/scratch/delete90/meshinchi_s/jlsmith3/GRCh37"
GTF="$TARGET/Reference_Data/GRCh37/gtf/Homo_sapiens.GRCh37.87.gtf"
GENOME="$TARGET/Reference_Data/GRCh37/fasta/genome/Homo_sapiens.GRCh37.dna.primary_assembly.fa"

#working directory
dir="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq/full_run1"
cd $dir


#Samples
fq=$1 # fq="cell1-hq_transcripts.collapsed.rep.fq"
samp=${fq%.fq}
expn=$2

# convert to fasta file - actually probably not necessary??
cat $fq | grep -E -A 1 "^@PB\.[0-9]" | grep -E -v "^--"| tr "@" ">" >${samp}.fasta

#create output directory
mkdir -p ${samp}

#Run sqanti_qc.py
sqanti_qc.py -x $GMAPDB -t 8 -o $samp -d $samp -fl ${samp%.rep}.abundance.txt ${samp}.fasta $GTF $GENOME 2> ${samp}_sqantiQC.stderr

#create output directory
mkdir -p ${samp}_filter

#run sqanti_filter.py
sqanti_filter.py -d ${samp}_filter -i $samp/${samp}_corrected.fasta $samp/${samp}_classification.txt 2> ${samp}_sqantiFilt.stderr

# usig a GTF formatted file cell1-hq_transcripts.collapsed.gff
# sqanti_qc.py -g -o ${samp}_GTF -d ${samp}_GTF \
# 		-fl ${samp%.rep}.abundance.txt ${samp%.rep}.gff $GTF $GENOME 2> ${samp}_sqantiQC_GTF.stderr &

#with Expression Data
# sqanti_qc.py -x $GMAPDB -t 8 -o $samp -d sqanti_results \
# 	-FL ${samp%.rep}.abundance.txt -e $expn \
# 	${samp}.fasta $GTF $GENOME
