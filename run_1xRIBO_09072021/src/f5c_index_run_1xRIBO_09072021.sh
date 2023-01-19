# activate conda environment
source /lustrehome/afonzino/anaconda3/bin/activate myenv

# change working directory
cd /lustre/bio_running/DNA_Ribo/run_1xRIBO_09072021/results

# launch f5c index command
f5c index -d /lustre/bio_running/DNA_Ribo/run_1xRIBO_09072021/workspace /lustre/bio_running/DNA_Ribo/run_1xRIBO_09072021/run_1xRIBO_09072021.fastq -t 10 --iop 10 --verbose 2> /lustre/bio_running/DNA_Ribo/run_1xRIBO_09072021/results/f5c_index_run_1xRIBO_09072021.err

# qsub command lauched
# qsub -q testqueue@sauron.recas.ba.infn.it -l nodes=1:ppn=15 f5c_index_run_1xRIBO_09072021.sh