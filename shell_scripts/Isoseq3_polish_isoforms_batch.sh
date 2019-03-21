#!/bin/bash

#Jenny Smith
#Run Isoseq3 polish

#EXAMPLE USAGE:
#BUCKET="s3://my-bucket"
#aws s3 cp isoseq3_polish_isoforms_batch.sh $BUCKET/my-scripts/
#aws batch submit-job file://submit_polish.json

#set script to exit 1 if any of the following are not met.
set -euo pipefail


#Define File Locations
BUCKET="s3://fh-pi-meshinchi-s"
SCRATCH="/scratch/${AWS_BATCH_JOB_ID}_${AWS_BATCH_CE_NAME}" #create unique directory on scratch volume

#Define Samples
PREFIX=${1:-"ccs_combined"} #if not given, use ccs_combined as the default PREFIX


echo $UNPOLISHED_BAM_PREFIX
echo $PREFIX
echo $AWS_BATCH_JOB_ID


#Copy the datasets into the container.
#This assumes all subreads.bam files in the specified BUCKET are included and were used to create the UNPOLISHED_BAM.
#The copying of the results of `isoseq3 cluster` to the container  should have no more than 3 files so those are specified.
echo "Copying the subreads.bam and unpolished.bam"
mkdir -p $SCRATCH
aws s3 cp --quiet --recursive --exclude "*" --include "*.subreads.bam"  $BUCKET/SR/SMRTSeq/ $SCRATCH
aws s3 cp --quiet --recursive --exclude "*" --include "${UNPOLISHED_BAM_PREFIX}.bam" --include "${UNPOLISHED_BAM_PREFIX}.bam.pbi" $BUCKET/SR/SMRTSeq/Isoseq3/ $SCRATCH

#Create a subread set for the original samples' "movie.subreads.bam" files.
#Use find command to create a list of the input subreads.bam files with full file paths.
#this assumes all subreads.bam files in the current working directory are included.
SUBREAD_BAMS=$(find $SCRATCH -type f -name "*.subreads.bam")
echo "Creating subreadset.xml"
dataset create --type SubreadSet --name $PREFIX ${PREFIX}.subreadset.xml $SUBREAD_BAMS

#run the serial polishing algorithm
echo "Running polish"
isoseq3 polish -j 32 --verbose --log-file "${PREFIX}_polish.log"  $SCRATCH/${UNPOLISHED_BAM_PREFIX}.bam ${PREFIX}.subreadset.xml ${PREFIX}.polished.bam

#Upload the results to the S3 bucket
aws s3 cp --quiet --recursive --exclude "*" --include "${PREFIX}*"  $PWD $BUCKET/SR/SMRTSeq/Isoseq3/

#remove the scratch directory
echo "Removing scratch directory."
rm -rf $SCRATCH
