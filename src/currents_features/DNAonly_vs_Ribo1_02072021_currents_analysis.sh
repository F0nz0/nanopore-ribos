# activate conda environment.
source /lustrehome/afonzino/anaconda3/bin/activate myenv

# define output and error paths
OUT_PATH=/lustre/bio_running/DNA_Ribo/src/DNAonly_vs_Ribo1_02072021_currents_analysis.sh.out
ERR_PATH=/lustre/bio_running/DNA_Ribo/src/DNAonly_vs_Ribo1_02072021_currents_analysis.sh.err

# lauch python script to produce graphs of currents
python3 /lustre/bio_running/DNA_Ribo/src/DNAonly_vs_Ribo1_02072021_currents_analysis.py >$OUT_PATH 2>$ERR_PATH

# qsub command
# qsub -q testqueue@sauron.recas.ba.infn.it -l nodes=1:ppn=20 DNAonly_vs_Ribo1_02072021_currents_analysis.sh