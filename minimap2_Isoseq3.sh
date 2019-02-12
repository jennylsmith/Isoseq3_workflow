#!/bin/bash
#SBATCH --partition=campus
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -o GMAP.%j.out
#SBATCH -e GMAP.%j.stderr
#SBATCH -t 0-1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlsmith3@fredhutch.org
source /app/Lmod/lmod/lmod/init/bash
################################

#Jenny Smith
#Minimap2 Alignment of Iso-seq3 Results for use with Cupcake Collapse Isoforms and count isoforms


DELETE90="/fh/scratch/delete90/meshinchi_s/jlsmith3"


############ STEP 1 #################
# Setup Ananconda
# https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda

# module purge
# ./Anaconda2-2018.12-Linux-x86_64.sh
# conda -V
# conda update conda

# conda create -n anaconda2.7 python=2.7 anaconda
# source activate anaconda2.7
# conda install -n  biopython
# conda install -n anaconda2.7 -c http://conda.anaconda.org/cgat bx-python
# conda install -n anaconda2.7 -c bioconda isoseq3
# conda install -n anaconda2.7 -c bioconda pbcoretools # for manipulating PacBio datasets
# conda install -n anaconda2.7 -c bioconda bamtools    # for converting BAM to fasta
# conda install -n anaconda2.7 -c bioconda pysam       # for making CSV


############## STEP 2 ################
#download precompiled binaries for minimap2
# curl -L https://github.com/lh3/minimap2/releases/download/v2.15/minimap2-2.15_x64-linux.tar.bz2 | tar -jxvf -



######## STEP 3 ###############
#for cDNA_cupcake with python2.7.6:
#https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#collapse
# git clone https://github.com/Magdoll/cDNA_Cupcake.git
# source activate anaconda2.7
# cd cDNA_Cupcake
# python setup.py build
# python setup.py install
# Installed /home/jlsmith3/anaconda2/envs/anaconda2.7/lib/python2.7/site-packages/bx_python-0.7.3-py2.7-linux-x86_64.egg


#Load Modules/define paths. Using python2.7 in virtual enviornment "anaconda2.7"
module purge
source activate anaconda2.7

export PATH=~/scripts/opt/bin:~/.local/bin:$PATH #for minimap2 & cDNA_cupcake
ml SAMtools/1.8-foss-2016b


#working directory
dir="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq"
GMAPDB="/fh/scratch/delete90/meshinchi_s/jlsmith3/GRCh37" #NEEDS TESTING
fq=$1 # fq="cell1-hq_transcripts.fastq"
samp=${fq%.fastq}


#Run alignment
gmap -D $GMAPDB -d GRCh37 -f samse -n 0 -t 6 --cross-species \
    --max-intronlength-ends 200000 -z sense_force $dir/$fq > $dir/${fq}.sam \
    2> ${samp}_gmap.stderr


#sort and index Bam/Sam files
sort -k 3,3 -k 4,4n  $dir/${fq}.sam > $dir/${fq}.srt.sam
samtools view -bS $dir/${fq}.sam > $dir/${fq}.bam
samtools sort $dir/${fq}.bam > $dir/${fq}.srt.bam
samtools index $dir/${fq}.srt.bam


#Collapse Reads
collapse_isoforms_by_sam.py --input $dir/${fq} --fq \
   -s $dir/${fq}.srt.sam --dun-merge-5-shorter -o $samp

#calculate abundance
get_abundance_post_collapse.py ${samp}.collapsed  $dir/${samp}-file.csv 2> ${samp}_counts.stderr

#Filter away 5' degraded isoforms
filter_away_subset.py ${samp}.collapsed



#GMAP Notes on Parameters.
   # --min-intronlength=INT  Min length for one internal intron (default 9).  Below this size,
       #                                    a genomic gap will be considered a deletion rather than an intron.
       #
   # -z, --direction=STRING  cDNA direction (sense_force, antisense_force,
       #                                                                     sense_filter, antisense_filter,or auto (default))
       # --cross-species     Use a more sensitive search for canonical splicing, which helps especially
       #
  # -n, --npaths=INT    Maximum number of paths to show (default 5).  If set to 1, GMAP
       #                                  will not report chimeric alignments, since those imply
       #                                  two paths.  If you want a single alignment plus chimeric
       #                                  alignments, then set this to be 0.
