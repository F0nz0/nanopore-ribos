# activate conda environment.
source /lustrehome/afonzino/anaconda3/bin/activate myenv

# launch minimap2 command to align bam to reference and to sort via samtools.
minimap2 -ax map-ont -t 8 /lustre/bio_running/DNA_Ribo/refs/ref.fa /lustre/bio_running/DNA_Ribo/run_dna_30062021/dna_only.fastq | /lustrehome/afonzino/samtools-1.13/samtools sort -o /lustre/bio_running/DNA_Ribo/run_dna_30062021/reads-ref_dna_only.sorted.bam -T /lustre/bio_running/DNA_Ribo/run_dna_30062021/reads.tmp 2> /lustre/bio_running/DNA_Ribo/run_dna_30062021/results/minimap2_dna_only.err

# launch samtools index command to index sorted bam file produced by the previous command.
/lustrehome/afonzino/samtools-1.13/samtools index /lustre/bio_running/DNA_Ribo/run_dna_30062021/reads-ref_dna_only.sorted.bam 2> /lustre/bio_running/DNA_Ribo/run_dna_30062021/results/samtools_index_dna_only.err

# check Let’s see if the bam file is not truncated. 
# This will check that the beginning of the file contains a valid header, and checks if the EOF is present. 
# This will exit with a non-zero exit code if the conditions were not met:
/lustrehome/afonzino/samtools-1.13/samtools quickcheck /lustre/bio_running/DNA_Ribo/run_dna_30062021/reads-ref_dna_only.sorted.bam> /lustre/bio_running/DNA_Ribo/run_dna_30062021/results/samtools_quickcheck_dna_only.out 2> /lustre/bio_running/DNA_Ribo/run_dna_30062021/results/samtools_quickcheck_dna_only.err

# qsub command lauched inside the /lustre/bio_running/DNA_Ribo/run_dna_30062021/src directory.
# qsub -q testqueue@sauron.recas.ba.infn.it -l nodes=1:ppn=20 minimap2_samtools_dna_only.sh