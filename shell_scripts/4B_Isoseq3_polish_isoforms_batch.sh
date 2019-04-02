#!/bin/bash

#Jenny Smith
#Run Isoseq3 polish

#EXAMPLE USAGE:
#BUCKET="s3://my-bucket"
#aws s3 cp 4A_isoseq3_polish_isoforms_batch.sh $BUCKET/my-scripts/
#aws batch submit-job file://submit_polish.json

#set script to exit 1 if any of the following are not met.
set -euo pipefail
export PATH=/home/ubuntu/anaconda2/bin:$PATH #re-declare path b/c fetch_and_run.sh overwrites PATH as the entry point.

#Define File Locations
BUCKET="s3://fh-pi-meshinchi-s"
SCRATCH="/scratch/${AWS_BATCH_JOB_ID}_${AWS_BATCH_CE_NAME}" #create unique directory on scratch volume
mkdir -p $SCRATCH
cd $SCRATCH #since `polish` will write output to the current working directory. So when it tries to do that in /home/ubuntu/ in the container it might crash/out of memory error.


#Define Samples based on the filename from `isoseq3 cluster`
PREFIX=$(echo ${UNPOLISHED_BAM_PREFIX}.bam | sed -E "s/\.bam|\.unpolished//gm")

echo $UNPOLISHED_BAM_PREFIX #this is the prefix used in 33_Isoseq3_cluster_isoforms.sh
echo $AWS_BATCH_JOB_ID

#Copy the datasets into the container.
#This assumes all subreads.bam files in the specified BUCKET are included and were used to create the UNPOLISHED_BAM.
echo "Copying the subreads.bam and unpolished.bam"
aws s3 cp --only-show-errors --recursive --exclude "*" --include "*.subreads.bam" --include "*.subreads.bam.pbi"  $BUCKET/SR/SMRTSeq/ .
aws s3 cp --only-show-errors --recursive --exclude "*" --include "${UNPOLISHED_BAM_PREFIX}.bam" --include "${UNPOLISHED_BAM_PREFIX}.bam.pbi" $BUCKET/SR/SMRTSeq/Isoseq3/ .
echo "$(ls -1R $PWD)"


#Create a subread set for the original samples' "movie.subreads.bam" files.
SUBREAD_BAMS=$(find $SCRATCH -type f -name "*.subreads.bam") #Use find command to create a list of the input subreads.bam files with full file paths.
echo "$SUBREAD_BAMS"
echo "Creating subreadset.xml"
dataset --debug create --force --type SubreadSet --name $PREFIX ${PREFIX}.subreadset.xml $SUBREAD_BAMS

#run the serial polishing algorithm
echo "Running polish"
isoseq3 polish --log-level "DEBUG" -j 32 --verbose   ${UNPOLISHED_BAM_PREFIX}.bam ${PREFIX}.subreadset.xml ${PREFIX}.polished.bam

#Upload the results to the S3 bucket
echo "Copying polish results to bucket"
aws s3 cp --only-show-errors --recursive --exclude "*" --include "${PREFIX}*" $PWD $BUCKET/SR/SMRTSeq/Isoseq3/

#remove the scratch directory
echo "Removing scratch directory"
cd /home/ubuntu/
rm -rf $SCRATCH
