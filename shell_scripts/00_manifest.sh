#!/bin/bash

#Jenny Smith
#February 11, 2019
#purpose: Create a manifest file of all samples sequenced
# and (maybe?) rename the samples to registration number for AML samples.


#Example Usage:


#set script to exit 1 if any of the following are not met.
set -euo pipefail

#Set-up environment
ml R/3.5.1-foss-2016b-fh1

#rawDir has the raw data file paths. It can be a single directory path or a list of directory file paths.
#rawDir="/fh/fast/meshinchi_s/pub/r54228_20190131_194324" or rawDir=$(find $PWD -type d -name "r5*")
rawDir="$*"

#function to get all full file paths of the subread.xml files
findpath () {
	find $1 -type f -name "*.subreadset.xml"
}

#define the output directory as the directory in which the first raw data set is saved.
OUTDIR="$(dirname $(echo $rawDir | cut -f 1 -d " "))"
export OUTDIR

#loop through the raw data directories.
for dir in $rawDir
do

  #loop through each subdirectory 1_A01, 2_B01, 3_C01, 4_D01 under the raw data directory.
  for file in $(findpath $dir)
  do
    reg=$(cat $file |  head -1 | cut -f  4 -d " " | sed -E 's/Name=.(.+)./\1/' )

    printf "$(dirname $file),$(basename $file),$reg\n" >> $OUTDIR/Sample_ID_Map.csv

    #rename each file with the patient registration number
    #for res in $(ls -1 $(dirname $file)/*)
    #do
    #  mv -n $res ${reg}_$res
    #  echo ${reg}_$res
    #done

  done

done

#Clean the sample manifest with TidyR
Rscript $(dirname $0)/0A_clean_sample_manifest.r

#remove the untidy version
rm $OUTDIR/Sample_ID_Map.csv
