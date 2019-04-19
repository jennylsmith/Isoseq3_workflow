#!/usr/bin/env python2


#Authors: Liz Tseng & Jenny Smith

#Purpose: Create a custom classify report (see `isoseq3 refine` for more info on how this was produced.)
#to define conditions for each movie pair, since a single sample condition was sequenced across 2 SMRT-cells (resulting in 2 movies).
#this classify report will  be used with demultiplexing step in cDNA cupcake ().

#import modules
import argparse
import os.path
import pandas as pd
from csv import DictReader, DictWriter

#set-up command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('flnc_report', help="Full File path to the flnc.report.csv file.")
parser.add_argument('manifest', help="Full File path to sample manifest produced with the 00_manifest.sh script.")
parser.add_argument('prefix', help="Prefix to append to the output modified flnc.report.csv file")
args = parser.parse_args()

#Sample Manifest to map movie subreads to sample condition
manifest=pd.read_csv(args.manifest)
manifest['Dataset'].replace(regex=True, to_replace=".subreadset.xml", value='', inplace=True)

#Dictionary to hold the movie ID to Sample condition mappings.
# d = {'m54228_181211_220100': 'NBM', 'm54228_181214_110428': 'NBM',
#      'm54228_190201_161538': 'AML', 'm54247_190125_201139': 'AML'}
d = manifest[['Dataset','Reg.']].set_index('Dataset').to_dict()['Reg.']

#create a file called flnc.report.hacked.csv for the new collapsed cluster report.
filename=os.path.dirname(args.flnc_report)+"/"+args.prefix+".flnc.report.hacked.csv"
f = open(filename, 'w')
reader = DictReader(open(args.flnc_report),delimiter=',') #open the original report.csv with each line as a dictionary
writer = DictWriter(f, reader.fieldnames, delimiter=',', lineterminator = '\n') # write to new file
writer.writeheader()

#loop through the flnc report and split by primer.
for r in reader:
    r['primer'] = d[r['id'].split('/')[0]]
    writer.writerow(r)

#close the hacked file
f.close()
