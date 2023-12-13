#!/bin/bash

INPUT_BAM=/lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/reads-ref_run_1xRIBO_02072021.reverse.sorted.bam
OUTPUT_BAM_BASEPATH=/lustre/bio_running/DNA_Ribo/mixed_datasets_whole_primer/RIBO
echo Input BAM: $INPUT_BAM

# 100K
echo "Performing 10k"
OUTPUT_BAM=${OUTPUT_BAM_BASEPATH}_100.reverse.bam
echo Output BAM: $OUTPUT_BAM
frac=$( samtools idxstats $INPUT_BAM | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=10000/total; if (frac > 1) {print 1} else {print frac}}' )
echo Fraction to extract: $frac
samtools view -bs $frac -r M13mp18:4975-5025 $INPUT_BAM > $OUTPUT_BAM
echo Computation Finished.

# 80K
echo "Performing 8k"
OUTPUT_BAM=${OUTPUT_BAM_BASEPATH}_80.reverse.bam
echo Output BAM: $OUTPUT_BAM
frac=$( samtools idxstats $INPUT_BAM | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=8000/total; if (frac > 1) {print 1} else {print frac}}' )
echo Fraction to extract: $frac
samtools view -bs $frac -r M13mp18:4975-5025 $INPUT_BAM > $OUTPUT_BAM
echo Computation Finished.

# 60K
echo "Performing 6k"
OUTPUT_BAM=${OUTPUT_BAM_BASEPATH}_60.reverse.bam
echo Output BAM: $OUTPUT_BAM
frac=$( samtools idxstats $INPUT_BAM | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=6000/total; if (frac > 1) {print 1} else {print frac}}' )
echo Fraction to extract: $frac
samtools view -bs $frac -r M13mp18:4975-5025 $INPUT_BAM > $OUTPUT_BAM
echo Computation Finished.

# 40K
echo "Performing 4k"
OUTPUT_BAM=${OUTPUT_BAM_BASEPATH}_40.reverse.bam
echo Output BAM: $OUTPUT_BAM
frac=$( samtools idxstats $INPUT_BAM | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=4000/total; if (frac > 1) {print 1} else {print frac}}' )
echo Fraction to extract: $frac
samtools view -bs $frac -r M13mp18:4975-5025 $INPUT_BAM > $OUTPUT_BAM
echo Computation Finished.

# 20K
echo "Performing 2k"
OUTPUT_BAM=${OUTPUT_BAM_BASEPATH}_20.reverse.bam
echo Output BAM: $OUTPUT_BAM
frac=$( samtools idxstats $INPUT_BAM | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=2000/total; if (frac > 1) {print 1} else {print frac}}' )
echo Fraction to extract: $frac
samtools view -bs $frac -r M13mp18:4975-5025 $INPUT_BAM > $OUTPUT_BAM
echo Computation Finished.

echo Total Computation Finished.