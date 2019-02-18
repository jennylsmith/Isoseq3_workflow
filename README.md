# Isoseq3_workflow


We have a total of 12 (11 AML, 1 Normal Bone Marrow) different biological samples to be sequenced across 26 SMRT-cells.  There are 4 SMRT-cells for the NBM.  Each AML sample is sequenced across 2 SMRT-cells (a large cDNA prep on 1 cell, and small cDNA prep on the 2nd cell).

We used 96 barcoded oligonucleotides for the RT-PCR step, to provide a means for PCR duplicate detection. Each SMRT-cell contains a single sample's library (not multiplexed by sample).

## General Workflow

### IsoSeq3
1. ccs generation
2. demultiplex the barcoded primers
3.

### Post-IsoSeq3 workflow
1. Align the *hq-transcripts.fasta files with GMAP-GSNAP
2. Sort/index the bams 
3. Use Cupcake ToFu to collapse the isoforms based on alignment coords
4. Use Cupcake ToFu to Calculate abundance of non-redundant isoforms
5. Run SQANTI2 using the unique transcripts from step #3 (*-hq_transcripts.collapsed.rep.fq )


## Resources and References

- https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki
- https://github.com/PacificBiosciences/IsoSeq_SA3nUP/wiki/Tutorial:-Installing-and-Running-Iso-Seq-3-using-Conda
- https://github.com/PacificBiosciences/IsoSeq3/blob/master/README_v3.1.md
- https://github.com/Magdoll/cDNA_Cupcake
- https://github.com/Magdoll/SQANTI2
- https://pacbiofileformats.readthedocs.io/en/3.0/BAM.html

