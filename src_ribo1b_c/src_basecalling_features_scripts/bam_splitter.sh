#!/bin/bash
# split forward and reverse reads generating two sorted and indexed bam files from the whole one.

BASEDIR=/lustre/bio_running/DNA_Ribo_2023
samtools=/lustrehome/afonzino/samtools-1.13/samtools
SAMPLE=$1
BAM=$BASEDIR/$SAMPLE/$SAMPLE.bam
BAM_REV=$BASEDIR/$SAMPLE/$SAMPLE.reverse.bam
BAM_FORW=$BASEDIR/$SAMPLE/$SAMPLE.forward.bam

echo Splitting BAM file in forward and reverse mapped reads:
echo Input BAM file: $BAM
echo Output Reverse BAM file: $BAM_REV
echo Ouput Forward BAM file: $BAM_FORW

echo Extracting Reverse reads and sorting these...
# reverse
# retrieve reverse reads
$samtools view -b -f 16 $BAM | $samtools sort -O BAM -o $BAM_REV
echo Indexing Reverse BAM file...
# index reverse bam file
$samtools index $BAM_REV

echo Extracting Forward reads and sorting these...
# forward
# retrieve forward reads
$samtools view -b -F 16 $BAM | $samtools sort -O BAM -o $BAM_FORW
echo Indexing Forward BAM file...
# index reverse bam file
$samtools index $BAM_FORW

echo Computation Finished.