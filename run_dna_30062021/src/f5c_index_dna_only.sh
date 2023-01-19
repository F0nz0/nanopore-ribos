# activate conda environment
source /lustrehome/afonzino/anaconda3/bin/activate myenv

# launch f5c index command
f5c index -d /lustre/bio_running/DNA_Ribo/run_dna_30062021/workspace /lustre/bio_running/DNA_Ribo/run_dna_30062021/dna_only.fastq -t 10 --iop 10 --verbose 2> /lustre/bio_running/DNA_Ribo/run_dna_30062021/results/f5c_index_dna_only.err

# qsub command lauched
# qsub -q testqueue@sauron.recas.ba.infn.it -l nodes=1:ppn=20 f5c_index_dna_only.sh