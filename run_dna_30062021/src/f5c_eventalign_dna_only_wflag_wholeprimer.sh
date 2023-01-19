# activate conda environment.
source /lustrehome/afonzino/anaconda3/bin/activate myenv

# f5c alignment with w flag on the whole region of the primer for dna only (the same of run_1xRIBO)
f5c eventalign --iop 9 -t 18 -r /lustre/bio_running/DNA_Ribo/run_dna_30062021/dna_only.fastq \
	-b /lustre/bio_running/DNA_Ribo/run_dna_30062021/reads-ref_dna_only.sorted.bam \
	-g /lustre/bio_running/DNA_Ribo/refs/ref.fa \
	--scale-events \
	--print-read-names \
	--samples \
	-w 'M13mp18:4970-5027' \
	--signal-index > /lustre/bio_running/DNA_Ribo/run_dna_30062021/f5c_eval_wflags/dna_only_wflag_wholeprimer.eventalign 2> /lustre/bio_running/DNA_Ribo/run_dna_30062021/results/f5c_eventalign_dna_only_wflag_wholeprimer.err

# qsub command lauched inside the /lustre/bio_running/DNA_Ribo/run_dna_30062021/src directory.
# qsub -q testqueue@sauron.recas.ba.infn.it -l nodes=1:ppn=20 f5c_eventalign_dna_only_wflag_wholeprimer.sh