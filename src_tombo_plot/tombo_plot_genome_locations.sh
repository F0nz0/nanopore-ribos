#!/bin/bash
# activate conda environment
source /lustrehome/afonzino/anaconda3/bin/activate tombo_env

# define inputs#
sleep_time=5

SRC_DIR=/lustre/bio_running/DNA_Ribo_minimap2sensitive/src_tombo_plot

OUTPUT_DIR=/lustre/bio_running/DNA_Ribo_minimap2sensitive/tombo_plot_outputs
mkdir $OUTPUT_DIR

SAMPLES=$1
FAST5SINGLE_DIR_CTRL=$2
FAST5SINGLE_DIR=$3
GENOME_LOCATION=$4
SITENAME=$5
STDERR=$SRC_DIR/tombo_plot_genome_locations.$SAMPLES.$SITENAME.err
PDF_FILEPATH_OUTPUT=$OUTPUT_DIR/tombo_results.$SAMPLES.$GENOME_LOCATION.$SITENAME.genome_location.pdf

echo "##############################################################################################" > $STDERR
echo Performing Tombo plot genome locations on samples $SAMPLES at site $SITENAME coordinates $GENOME_LOCATION >> $STDERR

# launch tombo plot genome_locations
tombo plot genome_locations \
        --fast5-basedirs $FAST5SINGLE_DIR \
        --control-fast5-basedirs $FAST5SINGLE_DIR_CTRL \
		--genome-locations $GENOME_LOCATION \
		--plot-standard-model \
		--overplot-threshold 50 \
		--overplot-type Density \
		--num-bases 21 \
		--pdf-filename $PDF_FILEPATH_OUTPUT 2>> $STDERR