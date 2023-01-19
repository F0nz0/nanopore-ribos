#!/bin/bash
# activate conda environment
source /lustrehome/afonzino/anaconda3/bin/activate tombo_env

# define inputs
sleep_time=5

SRC_DIR=/lustre/bio_running/DNA_Ribo_minimap2sensitive/src_tombo

SAMPLE=$1
FAST5SINGLE_DIR=$2
REF=/lustre/bio_running/DNA_Ribo/refs/ref.fa
STDERR=$SRC_DIR/tombo_resquiggle.$SAMPLE.err

echo "##############################################################################################" > $STDERR
echo Performing Tombo resquiggle on sample $SAMPLE >> $STDERR

tombo resquiggle \
	$FAST5SINGLE_DIR \
	$REF \
	--dna \
	--basecall-group Basecall_1D_005 \
	--overwrite \
	--ignore-read-locks \
	--failed-reads-filename $SRC_DIR/tombo_resquiggle_failed_reads.$SAMPLE \
	--processes 60 2>> $STDERR