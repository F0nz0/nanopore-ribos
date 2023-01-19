# activate conda environment.
source /lustrehome/afonzino/anaconda3/bin/activate myenv

# change working directory
cd /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/results

# launch minimap2 command to align bam to reference and to sort via samtools.
minimap2 -ax map-ont -t 8 /lustre/bio_running/DNA_Ribo/refs/ref.fa /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/run_1xRIBO_02072021.fastq | /lustrehome/afonzino/samtools-1.13/samtools sort -o /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/reads-ref_run_1xRIBO_02072021.sorted.bam -T /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/reads.tmp 2> /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/results/minimap2_run_1xRIBO_02072021.err

# launch samtools index command to index sorted bam file produced by the previous command.
/lustrehome/afonzino/samtools-1.13/samtools index /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/reads-ref_run_1xRIBO_02072021.sorted.bam 2> /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/results/samtools_index_run_1xRIBO_02072021.err

# check Let’s see if the bam file is not truncated. 
# This will check that the beginning of the file contains a valid header, and checks if the EOF is present. 
# This will exit with a non-zero exit code if the conditions were not met:
/lustrehome/afonzino/samtools-1.13/samtools quickcheck /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/reads-ref_run_1xRIBO_02072021.sorted.bam> /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/results/samtools_quickcheck_run_1xRIBO_02072021.out 2> /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/results/samtools_quickcheck_run_1xRIBO_02072021.err

# qsub command lauched inside the /lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/src directory.
# qsub -q testqueue@sauron.recas.ba.infn.it -l nodes=1:ppn=10 minimap2_samtools_run_1xRIBO_02072021.sh