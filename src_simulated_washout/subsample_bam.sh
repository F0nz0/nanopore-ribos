#!/bin/bash

INPUT_BAM=$1
OUTPUT_BAM=$2
n_reads=$3
echo Input BAM: $INPUT_BAM
echo Output BAM: $OUTPUT_BAM
echo Number of reads to subsample: $n_reads
frac=$( samtools idxstats $INPUT_BAM | cut -f3 | awk -v N_READS="$n_reads" 'BEGIN {total=0} {total += $1} END {frac=N_READS/total; if (frac > 1) {print 1} else {print frac}}' )
echo Fraction to extract: $frac
samtools view -bs $frac $INPUT_BAM > $OUTPUT_BAM
echo Computation Finished.