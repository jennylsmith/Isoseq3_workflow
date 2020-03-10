#!/bin/bash


#Jenny Smith
#Feb 13, 2019
#purpose: create MD5 sums recursively.


#initial directory. this is a positional argument - must type the filepath on the command line.
dir1="$1"

#move to original files directory
cd $dir1


findpath () {
	find $1 -type f -regex "^.+subreads\.bam$"
}

#Loop through each directory with raw SMRT-seq data
for dir in $(ls -1d r5* )
do
	subreads_bam=$( findpath $dir )

	#create md5 sums for each subreads bam.
	for bam in $(echo $subreads_bam)
	do
		echo $bam
		lastDir=$(dirname $bam | uniq)
		sample=$(basename $bam)
		m=$(cd $lastDir; md5sum $sample > ${sample}.MD5)
	done

done
