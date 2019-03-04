#!bin/bash

#Jenny Smith
#February 11, 2019
#purpose: rename the samples to registration number for AML samples.


#rawDir="/fh/fast/meshinchi_s/pub/r54228_20190131_194324"
rawDir=$1
mappingFile=$2

findpath () {
	find $1 -type f -name "*.subreadset.xml"
}

#loop through each subdirectory 1_A01, 2_B01, 3_C01, 4_D01 under the raw data directory.
for file in $(findpath $rawDir)
do
 reg=$(cat $file |  head -1 | cut -f  4 -d " " | sed -E 's/Name=.(.+)./\1/' )

 #rename each file with the patient registration number
 #for res in $(ls -1 $(dirname $file)/*)
 #do
 # mv -n $res ${reg}_$res
  #echo ${reg}_$res
 #done

 printf "$(dirname $file),$(basename $file),$reg\n"
done #>> Sample_ID_Map.csv
