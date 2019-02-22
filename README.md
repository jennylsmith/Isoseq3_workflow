# Isoseq3_workflow


We have a total of 12 (11 AML, 1 Normal Bone Marrow) different biological samples to be sequenced across 26 SMRT-cells.  There are 4 SMRT-cells for the NBM.  Each AML sample is sequenced across 2 SMRT-cells (a large cDNA prep on 1 cell, and small cDNA prep on the 2nd cell).

We used 96 barcoded oligonucleotides for the RT-PCR step, to provide a means for PCR duplicate detection. Each SMRT-cell contains a single sample's library (not multiplexed by sample).

## General Workflow

### IsoSeq3
1. ccs generation
2. demultiplex the barcoded primers
3.

### Post-IsoSeq3 workflow
1. Align the *hq-transcripts.fastq files with minimap2
2. Sort/index the aligned BAM files 
3. Use Cupcake ToFu to collapse the isoforms based on alignment coords
4. Use Cupcake ToFu to Calculate abundance of non-redundant isoforms
5. Run SQANTI2 using the unique transcripts from step #3 (*-hq_transcripts.collapsed.rep.fq)


### Example Workflow (Liz)

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
isoseq3 refine Jenny_4cell_FL.consensusreadset.xml barcodes.fasta flnc.bam
```
(to be continued)


## Resources and References

- https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki
- https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda
- https://github.com/PacificBiosciences/IsoSeq3/blob/master/README_v3.1.md
- https://github.com/Magdoll/cDNA_Cupcake
- https://github.com/Magdoll/SQANTI2
- https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
- https://github.com/PacificBiosciences/PacBioFileFormats/wiki/BAM-recipes 

