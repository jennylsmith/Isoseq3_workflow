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
#GMAP-SNAP Alignment of Iso-seq3 Results for use with Cupcake Collapse Isoforms and count isoforms


DELETE90="/fh/scratch/delete90/meshinchi_s/jlsmith3"

######## STEP 1 ###############
#Version of GMAP-GSNAP to install
#wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-07-04.tar.gz

#Edited the config file
#config.site: prefix=/home/jlsmith3/scripts/opt, with_gmapdb=/fh/scratch/delete90/meshinchi_s/GRCh37, PERL=/app/easybuild/software/Perl/5.28.0-foss-2016b

#Used below to install gmap-2018-07-04/
# ./configure
# make
# make check
# make install
# gmap_build -D $DELETE90 -d GRCh37 -g $(cat fasta_files.txt) 2> gmap_build2.stderr &
# make gmapdb  #NOTE:had a weird error. states "make: no recipe for gmapdb" ??? but appears to be working



######## STEP 2 ###############
#for cDNA_cupcake with python2.7.6:
#https://github.com/Magdoll/cDNA_Cupcake/wiki/Cupcake-ToFU:-supporting-scripts-for-Iso-Seq-after-clustering-step#collapse
# pip install --user bx-python
# pip install --user biopython
# git clone https://github.com/Magdoll/cDNA_Cupcake.git
# cd cDNA_Cupcake
# python setup.py build
# python setup.py install --user
# python -c "import site; print(site.USER_BASE)" #(where my user scheme installs programs)


#Load Modules/define paths. Python must be /usr/bin/python for python 2.7.6
module purge
export PATH=~/scripts/opt/bin:~/.local/bin:$PATH #for GMAP & cDNA_cupcake
ml SAMtools/1.8-foss-2016b


#working directory
dir="/fh/scratch/delete90/meshinchi_s/jlsmith3/SMRTseq"
GMAPDB="/fh/scratch/delete90/meshinchi_s/jlsmith3/GRCh37" #NEEDS TESTING
fq=$1 # fq="cell1-hq_transcripts.fastq"
flnc=$2  #fasta=cell1-unpolished.flnc.bam
samp=${fq%.fastq}

#Convert FLNC bam to fasta??
# flnc_id=${flnc%.bam}
# samtools fasta -@ 8 $flnc > ${flnc%.bam}.fasta &

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


#Next steps???
# Filter collapse results by minimum FL count support
# filter_by_count.py test.collapsed --min_count=2

#Filter away 5' degraded isoforms
# filter_away_subset.py test.collapsed.min_fl_2



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
