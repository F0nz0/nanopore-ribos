#!/bin/bash
# activate conda environment
source /lustrehome/afonzino/anaconda3/bin/activate tombo_env

# define inputs
sleep_time=5

SRC_DIR=/lustre/bio_running/DNA_Ribo_minimap2sensitive/src_tombo

SAMPLE=$1
FAST5DIR=$2
SAVEDIR=$3
STDERR=$SRC_DIR/$SAMPLE.err

echo "##############################################################################################" > $STDERR
echo Performing multi_to_single_fast5 program on sample $SAMPLE >> $STDERR

multi_to_single_fast5 --input_path $FAST5DIR --save_path $SAVEDIR -t 5 --recursive 2>> $STDERR