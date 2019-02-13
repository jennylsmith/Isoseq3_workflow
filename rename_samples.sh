#!bin/bash

#Jenny Smith
#February 11, 2019
#purpose: rename the samples to registration number for AML samples.



#mappingFile="/Users/work/Documents/GitHub/Isoseq3_workflow/samples/ibraries_submitted_to_FHCRC_genomics_core_for_PacBio.csv"
#rawDir="r54228_20190131_194324"
rawDir=$1
mappingFile=$2

#loop through each subdirectory 1_A01, 2_B01, 3_C01, 4_D01 under the raw data directory.
for file in $(ls -1 $rawDir/{1..4}_*01/*.subreadset.xml)
do
 reg=$(cat $file |  head -1 | cut -f  4 -d " " | sed -E 's/Name=.(.+)./\1/' )

 #rename each file with the patient registration number
 #for res in $(ls -1 $(dirname $file)/*)
 #do
 # mv -n $res ${reg}_$res
  #echo ${reg}_$res
 #done

 printf "$(dirname $file),$(basename $file),$reg\n"
done >> Sample_ID_Map.txt
