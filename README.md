# Isoseq3_workflow

This is a custom PacBio RNA sequencing workflow, which is a full-length RNA-seq methodology that can sequence transcripts ranging from 10 kb and up. 

We have a total of 11 (10 AML, 1 Normal Bone Marrow) biological samples sequenced across 26 SMRT-cells.  There are 4 SMRT-cells for the noraml bone marrow (NBM) sample. Each AML sample is sequenced across 2 SMRT-cells (a large cDNA prep on 1 cell, and small cDNA prep on the 2nd cell).

We used 96 barcoded oligonucleotides for the RT-PCR step during library preparation, to provide a means for PCR duplicate detection. Each SMRT-cell contains a single sample's library (not multiplexed by sample).


## General Workflow

This workflow has been initially designed to be run using the SLURM job scheduler with an on-premise high performance compute cluster (HPC). 

The next iteration will use Cromwell Workflow Management System to improve the efficiency and reprodicibility of the workflow. 

### IsoSeq3

1. CCS (Circular Consensus Sequence) generation for each individual movie. 

```
for bam in $(ls -1 *.subreads.bam) 
do
  sbatch Isoseq3_ccs_generation.sh $bam
done
```

2. Create a ConsensusReadSet of each movie's ccs.bam file. Steps 2-5 are completed in a single batch script. 
   - This creates an xml file which provides the instructions on which datasets to include in a particular command. 
   - This will allow one to combine multiple movies in the `lima` command  in order to remove 5' and 3' primers and demultiplex barcoded samples. 

3. Use the `lima` command on the combined .ccs.bam files.
   - The .ccs.bams are listed in the .xml file from step 2 and specifiy to `lima` the location of the .ccs.bam files. The actual merging/combining of the .ccs.bam files are handled by `lima`.  

4. Create a ConsensusReadSet from the output of `lima` which will help improve transcript recovery in the later clustering steps. 

5. Use `isoseq3 refine` to remove polyA tails and artificial concatemers from the combined demultiplexed bams produced by `lima`. 
   - The refine step produces the full-length, non-concatemer (FLNC) reads. 
   - The bam files can be listed from the command line - in this workflow, I assume all .ccs.bams in the currenty working directory will be used. 
   - You will need the barcode and primer sequences used in the library prep in fasta format for this step. In this workflow, the barcodes fasta file is aptly named `barcodes.fasta`. 
   - The prefix is simply a string to identify your run like "test1", or "AML". 

```
ccs_bams=$(ls *.ccs.bam)
prefix="Test"
sbatch Isoseq3_lima_primerRemoval.sh $ccs_bams barcodes.fasta $prefix 

```

6. Clustering of the FLNC reads in `isoseq3 cluster` 

   - clustering is the initial isoform-level clustering of full-length reads. This is the step which requires more power for isoform detection and thus, why combining as many samples in the initial steps 2-5 were necessary. 
   - Iso-Seq 3 requires at least **two FLNC reads** to b clustered at the isoform level
   - The two reads must *A)* Differ less than 100 bp on the 5’ end, *B)*  Differ less than 30 bp on the 3’ end, and *c)* Have no internal gaps that exceed 10 bp to be clustered. 


### Post-IsoSeq3 workflow

1. Align the *hq-transcripts.fastq files with minimap2
2. Sort/index the aligned BAM files 
3. Use Cupcake ToFu to collapse the isoforms based on alignment coords
4. Use Cupcake ToFu to Calculate abundance of non-redundant isoforms
5. Run SQANTI2 using the unique transcripts from step #3 (*-hq_transcripts.collapsed.rep.fq)


### Example Workflow (Liz)

Last Updated: 03/25/2019

Liz uses 4 cells to demonstrate how to first run LIMA and IsoSeq3 on combined cells (for higher throughput) and later demultiplex into per-sample isoform counts.

1. Combine the CCS into a `consensusreadset` and run LIMA

```
# lima 1.8.0 (commit v1.8.0)
dataset create --type ConsensusReadSet --name Jenny_4cell Jenny_4cell.consensusreadset.xml m54228_181211_220100.ccs.bam m54228_181214_110428.ccs.bam m54228_190201_161538.ccs.bam m54247_190125_201139.ccs.bam
lima --isoseq --dump-clips Jenny_4cell.consensusreadset.xml barcodes.fasta output.bam
```

Validation `output.lima.summary`
```
ZMWs input                (A) : 1930352
ZMWs above all thresholds (B) : 1412342 (73%)
ZMWs below any threshold  (C) : 518010 (27%)
```

2. Combine the LIMA demux output into a `consensusreadset`

```
ls output.primer_5p--BC0*bam|xargs echo "dataset create --type ConsensusReadSet --name Jenny_4cell_FL Jenny_4cell_FL.consensusreadset.xml  " > cmd
bash cmd
```

3. Run IsoSeq3

```
isoseq3 refine Jenny_4cell_FL.consensusreadset.xml barcodes.fasta flnc.bam --require-polya
isoseq3 cluster -j 16 --verbose --log-file "test_cluster.log"  test.flnc.bam test.unpolished.bam --split-bam 24

# can run below commands in parallel
isoseq3 polish unpolished.0.bam Jenny_4cell.subreadset.xml polished.0.bam -j 12
isoseq3 polish unpolished.1.bam Jenny_4cell.subreadset.xml polished.1.bam -j 12
...
isoseq3 polish unpolished.23.bam Jenny_4cell.subreadset.xml polished.23.bam -j 12

# then combine the HQ FASTQ
gunzip polished.*.hq.fastq.gz
cat polished.18.hq.fastq > test_polished.hq.fastq
```

4. Map to Genome and Collapse

See [Cupcake tutorial here](https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/What-to-do-after-Iso-Seq-Cluster%3F)

```
minimap2 -ax splice -t 30 -uf --secondary=no hg38_noalt.fa test_polished.hq.fasta > test_polished.hq.fasta.sam
sort -k 3,3 -k 4,4n test_polished.hq.fasta.sam > test_polished.hq.fasta.sorted.sam
collapse_isoforms_by_sam.py --input test_polished.hq.fastq --fq -s test_polished.hq.fasta.sorted.sam -c 0.99 -i 0.95 --dun-merge-5-shorter -o test_polished.hq.no5merge

get_abundance_post_collapse.py test_polished.hq.no5merge.collapsed test_polished.cluster_report.csv
filter_away_subset.py test_polished.hq.no5merge.collapsed
```

5. Hack classify report and demux into two conditions

Run the following Python script to hack by movie to change the `primer` field to reflect the two conditions.


```
from csv import DictReader, DictWriter

d = {'m54228_181211_220100': 'NBM', 'm54228_181214_110428': 'NBM', \
     'm54228_190201_161538': 'AML', 'm54247_190125_201139': 'AML'}
     
f = open('flnc.report.hacked.csv', 'w')
reader = DictReader(open('flnc.report.csv'),delimiter=',')
writer = DictWriter(f, reader.fieldnames, delimiter=',')
writer.writeheader()
for r in reader:
    r['primer'] = d[r['id'].split('/')[0]]
    writer.writerow(r)

f.close()
```

Next run the demux script. See [tutorial here](https://github.com/Magdoll/cDNA_Cupcake/wiki/Tutorial:-Demultiplexing-SMRT-Link-Iso-Seq-Jobs#mapped)

```
python <path_to_cupcake>/post_isoseq_cluster/demux_isoseq_with_genome.py --mapped_fastq test_polished.hq.no5merge.collapsed.filtered.rep.fq --read_stat test_polished.hq.no5merge.collapsed.read_stat.txt --classify_csv ../flnc.report.hacked.csv -o test_polished.hq.no5merge.collapsed.filtered.mapped_fl_count.txt
```

## Resources and References

- https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki
- https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda
- https://github.com/PacificBiosciences/IsoSeq3/blob/master/README_v3.1.md
- https://github.com/Magdoll/cDNA_Cupcake
- https://github.com/Magdoll/SQANTI2
- https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
- https://github.com/PacificBiosciences/PacBioFileFormats/wiki/BAM-recipes 
- https://pacbiofileformats.readthedocs.io/en/5.1/


