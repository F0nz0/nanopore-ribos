#!/bin/bash

INPUT_BASEPATH_DNA_BAM=/lustre/bio_running/DNA_Ribo/mixed_datasets_whole_primer/DNA
INPUT_BASEPATH_RIBO_BAM=/lustre/bio_running/DNA_Ribo/mixed_datasets_whole_primer/RIBO
RIBO_NAME=ribo1
OUTPUT_BASEPATH_BAM=/lustre/bio_running/DNA_Ribo/mixed_datasets_whole_primer

echo Input DNA BAM BASEPATH: $INPUT_BASEPATH_DNA_BAM
echo Input RIBO BAM BASEPATH: $INPUT_BASEPATH_RIBO_BAM

# Performing mergings
for perc in 80 60 40 20;
    do
        ribo_perc=$perc
        dna_perc=$(expr 100 - $ribo_perc)
        echo ""
        echo Performing $dna_perc:$ribo_perc - DNA:RIBO
        INPUT_DNA_BAM=${INPUT_BASEPATH_DNA_BAM}_${dna_perc}.reverse.bam
        INPUT_RIBO_BAM=${INPUT_BASEPATH_RIBO_BAM}_${ribo_perc}.reverse.bam
        OUTPUT_BAM=$OUTPUT_BASEPATH_BAM/mix_dna_${RIBO_NAME}_${dna_perc}_${ribo_perc}.reverse.bam
        echo Input DNA BAM: $INPUT_DNA_BAM
        echo Input RIBO BAM: $INPUT_RIBO_BAM
        echo Output BAM: $OUTPUT_BAM
        samtools merge -O BAM $OUTPUT_BAM $INPUT_DNA_BAM $INPUT_RIBO_BAM
        echo Computation Finished for current mix.
    done

echo ""
echo Producing mix_dna_ribo_0_100
INPUT_BAM=${INPUT_BASEPATH_RIBO_BAM}_100.reverse.bam
OUTPUT_BAM=$OUTPUT_BASEPATH_BAM/mix_dna_${RIBO_NAME}_0_100.reverse.bam
echo Input BAM: $INPUT_BAM
echo Output file: $OUTPUT_BAM
cp $INPUT_BAM $OUTPUT_BAM

echo ""
echo Producing mix_dna_ribo_100_0
INPUT_BAM=${INPUT_BASEPATH_DNA_BAM}_100.reverse.bam
OUTPUT_BAM=$OUTPUT_BASEPATH_BAM/mix_dna_${RIBO_NAME}_100_0.reverse.bam
echo Input BAM: $INPUT_BAM
echo Output file: $OUTPUT_BAM
cp $INPUT_BAM $OUTPUT_BAM

echo ""
echo Performing indexing of mixed merged BAM files into directory: $OUTPUT_BASEPATH_BAM
for bam in $(ls $OUTPUT_BASEPATH_BAM/mix_dna_${RIBO_NAME}*bam);
    do
        echo Indexing $bam
        samtools index $bam
    done

echo ""
echo Computation Finished! Exiting...