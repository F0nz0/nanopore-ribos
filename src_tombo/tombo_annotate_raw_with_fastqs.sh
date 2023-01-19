#!/bin/bash
# activate conda environment
source /lustrehome/afonzino/anaconda3/bin/activate tombo_env

# define inputs
sleep_time=5

SRC_DIR=/lustre/bio_running/DNA_Ribo_minimap2sensitive/src_tombo

SAMPLE=$1
FAST5SINGLE_DIR=$2
FASTQFILE=$3
STDERR=$SRC_DIR/tombo_annotate_raw_with_fastqs.$SAMPLE.err

echo "##############################################################################################" > $STDERR
echo Annotating fast5 single files with guppy fastq basecalling files $SAMPLE >> $STDERR

tombo preprocess annotate_raw_with_fastqs \
	--fast5-basedir $FAST5SINGLE_DIR \
	--fastq-filenames $FASTQFILE \
	--basecall-group Basecall_1D_005 \
	--processes 20 2>> $STDERR