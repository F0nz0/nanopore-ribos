# for RIBO1 rep1 02/07/2021

# importing basic modules
import pandas as pd
import pysam
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn
import os

# to plot better figures
plt.tight_layout()

# ---------- defining Utils Functions --------------- #

def get_results_list(bam_file_path, contig, start_pos_1_based, end_pos_1_based):
    bam_file = pysam.AlignmentFile(bam_file_path)
    start_pos_0_based = start_pos_1_based - 1
    end_pos_0_based = end_pos_1_based
    
    result = []
    for pileupcolumn in bam_file.pileup(contig, start_pos_0_based, end_pos_0_based, truncate=True, min_base_quality=0, max_depth=10000000):
        result.append(pileupcolumn.get_query_sequences(mark_matches=True, add_indels=True))
    bam_file.close()
    
    return result


def get_stats(result_list_single_pos, verbose=True):
    T_count = 0
    C_count = 0
    G_count = 0
    A_count = 0
    del_count = 0
    for i in result_list_single_pos:
        if i.lower() == "t":
            T_count += 1
        elif i.lower() == "c":
            C_count += 1
        elif i.lower() == "g":
            G_count += 1
        elif i.lower() == "a":
            A_count += 1
        if i.lower() == "*":
            del_count += 1
        if "+" in i.lower():
            i_splitted = i.lower().split("+")
            if i_splitted[0] == "t":
                T_count += 1
            elif i_splitted[0] == "c":
                C_count += 1
            elif i_splitted[0] == "g":
                G_count += 1
            elif i_splitted[0]  == "a":
                A_count += 1
        if "-" in i.lower():
            i_splitted = i.lower().split("-")
            if i_splitted[0] == "t":
                T_count += 1
            elif i_splitted[0] == "c":
                C_count += 1
            elif i_splitted[0] == "g":
                G_count += 1
            elif i_splitted[0]  == "a":
                A_count += 1
                
    if verbose:
        print("T", T_count)
        print("C", C_count)
        print("G", G_count)
        print("A", A_count)
        print("Deletions", del_count)
    
    output = [T_count, C_count, G_count, A_count, del_count]
    return output


def get_insertions(bam_file_path, contig, start_pos_1_based, end_pos_1_based):
    bam_file = pysam.AlignmentFile(bam_file_path)
    start_pos_0_based = start_pos_1_based - 2 # -2 cause we want to calculate insertion of the previous reference positions in 0-based index
    end_pos_0_based = end_pos_1_based - 1 # -1 cause we want to calculate insertion of the previous reference positions in 0-based index with right opened intervals (pythonic way)
    
    start_time = datetime.now()
    result = []
    for pileupcolumn in bam_file.pileup(contig, start_pos_0_based, end_pos_0_based, truncate=True, min_base_quality=0, max_depth=10000000):
        result.append(pileupcolumn.get_query_sequences(mark_matches=True, add_indels=True))
    stop_time = datetime.now()
    bam_file.close()
    
    final_output_insertions = []
    for pos in result:
        ins_count = 0
        for read in pos:
            if "+" in read.lower():
                ins_count += 1
        final_output_insertions.append(ins_count)
    return final_output_insertions


def calc_freq_bases(df):
    freq = pd.DataFrame(df.iloc[:,0:6].values / np.array(df.iloc[:,0:6].sum(axis=1).values).reshape(df.shape[0],1), 
             columns=["T_freq", "C_freq", "G_freq", "A_freq", "del_freq", "ins_freq"])
    freq["index"] = df.index.to_list()
    freq.set_index("index", inplace=True)
    df_with_freq = pd.concat([df, freq], axis=1)
    return df_with_freq


def get_stats_region(result_list_region, start_pos_1_based=None, stop_pos_1_based=None):
    final_output = []
    for i in result_list_region:
        output = get_stats(i, verbose=False)
        final_output.append(output)
    df_final_output = pd.DataFrame(final_output, columns=["T_count", "C_count", "G_count", "A_count", "del_count"])
    if start_pos_1_based != None and stop_pos_1_based != None:
        index = range(start_pos_1_based, stop_pos_1_based+1)
        df_final_output.index = index
    return df_final_output


def get_stats_from_region(bam_file_path, contig, start_pos_1_based, end_pos_1_based):
    start_time = datetime.now()
    r = get_results_list(bam_file_path, contig, start_pos_1_based, end_pos_1_based)
    df_r = get_stats_region(r, start_pos_1_based, end_pos_1_based)
    insertions = get_insertions(bam_file_path, contig, start_pos_1_based, end_pos_1_based)
    df_r["ins_freq"] = insertions
    df_r_freq = calc_freq_bases(df_r)
    stop_time = datetime.now()
    print(f"[{datetime.now()}] Elapsed time: {stop_time - start_time}", flush=True)
    return df_r_freq


# defining input files paths
ref_path = "/lustre/bio_running/DNA_Ribo/refs/ref.fa"
dna_rev_path = "/lustre/bio_running/DNA_Ribo/run_dna_30062021/reads-ref_dna_only.reverse.sorted.bam"
ribo_rev_path = "/lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/reads-ref_run_1xRIBO_02072021.reverse.sorted.bam"
dna_forw_path =  "/lustre/bio_running/DNA_Ribo/run_dna_30062021/reads-ref_dna_only.forward.sorted.bam"
ribo_forw_path = "/lustre/bio_running/DNA_Ribo/run_1xRIBO_02072021/reads-ref_run_1xRIBO_02072021.forward.sorted.bam"
results_folder_path = "/lustre/bio_running/DNA_Ribo/results/basecalling_features_with_raw_data"

# -------------- Start Computations ----------------#
start_time = datetime.now()
print(f"[{datetime.now()}] Starting computations.", flush=True)

# SITE 1 (rA)
site = 4985
start = site - 7
stop = site + 7
suptitle = f"DNA only vs Ribo1 site M13mp18:{site}"
print(f"[{datetime.now()}] Processing site {suptitle}.", flush=True)

# dna only forward
df_dna_forw_s1 = get_stats_from_region(dna_forw_path, "M13mp18", start, stop)
# dna only reverse
df_dna_rev_s1 = get_stats_from_region(dna_rev_path, "M13mp18", start, stop)
# Ribo forward
df_ribo_forw_s1 = get_stats_from_region(ribo_forw_path, "M13mp18", start, stop)
# Ribo reverse
df_ribo_rev_s1 = get_stats_from_region(ribo_rev_path, "M13mp18", start, stop)

# PLOTTING RESULTS
f, axes = plt.subplots(nrows=2, ncols=2, figsize=(15,7))
# plot results of dna only forward strand
df_dna_forw_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[0,0], legend=None,
                                color=["red", "blue", "orange", "green", "gray", "black"])
axes[0,0].set_title("DNA_only (+) strand")
axes[0,0].set_xticks([])
# plot results of ribo forward strand
df_ribo_forw_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[0,1], legend=None,
                                  color=["red", "blue", "orange", "green", "gray", "black"])
axes[0,1].set_title("Ribo (+) strand")
axes[0,1].set_xticks([])
# plot results for dna only reverse strand
df_dna_rev_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[1,0], legend=None,
                               color=["red", "blue", "orange", "green", "gray", "black"])
axes[1,0].set_title("DNA_only (-) strand")
# plot results for Ribo reverse strand
df_ribo_rev_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[1,1], legend=None, 
                                 color=["red", "blue", "orange", "green", "gray", "black"])
axes[1,1].set_title("Ribo (-) strand")
f.suptitle(suptitle)
axes[0,0].set_ylabel("Frequency")
axes[1,0].set_ylabel("Frequency")
axes[1,0].set_xlabel("Reference Position")
axes[1,1].set_xlabel("Reference Position")
f.legend(["T_freq", "C_freq", "G_freq", "A_freq", "Del_freq", "Ins_freq"], loc="right")
f.savefig(os.path.join(results_folder_path, f"{suptitle}.tiff"), dpi=300)

# save to disk
df_dna_forw_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.dna.forw.tsv"), sep="\t")
df_dna_rev_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.dna.rev.tsv"), sep="\t")
df_ribo_forw_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.ribo.forw.tsv"), sep="\t")
df_ribo_rev_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.ribo.rev.tsv"), sep="\t")

# ----- Start Computations ----------------#
# SITE 2 (rG)
site = 4997
start = site - 7
stop = site + 7
suptitle = f"DNA only vs Ribo1 site M13mp18:{site}"
print(f"[{datetime.now()}] Processing site {suptitle}.", flush=True)

# dna only forward
df_dna_forw_s1 = get_stats_from_region(dna_forw_path, "M13mp18", start, stop)
# dna only reverse
df_dna_rev_s1 = get_stats_from_region(dna_rev_path, "M13mp18", start, stop)
# Ribo forward
df_ribo_forw_s1 = get_stats_from_region(ribo_forw_path, "M13mp18", start, stop)
# Ribo reverse
df_ribo_rev_s1 = get_stats_from_region(ribo_rev_path, "M13mp18", start, stop)

# PLOTTING RESULTS
f, axes = plt.subplots(nrows=2, ncols=2, figsize=(15,7))
# plot results of dna only forward strand
df_dna_forw_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[0,0], legend=None,
                                color=["red", "blue", "orange", "green", "gray", "black"])
axes[0,0].set_title("DNA_only (+) strand")
axes[0,0].set_xticks([])
# plot results of ribo forward strand
df_ribo_forw_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[0,1], legend=None,
                                  color=["red", "blue", "orange", "green", "gray", "black"])
axes[0,1].set_title("Ribo (+) strand")
axes[0,1].set_xticks([])
# plot results for dna only reverse strand
df_dna_rev_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[1,0], legend=None,
                               color=["red", "blue", "orange", "green", "gray", "black"])
axes[1,0].set_title("DNA_only (-) strand")
# plot results for Ribo reverse strand
df_ribo_rev_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[1,1], legend=None, 
                                 color=["red", "blue", "orange", "green", "gray", "black"])
axes[1,1].set_title("Ribo (-) strand")
f.suptitle(suptitle)
axes[0,0].set_ylabel("Frequency")
axes[1,0].set_ylabel("Frequency")
axes[1,0].set_xlabel("Reference Position")
axes[1,1].set_xlabel("Reference Position")
f.legend(["T_freq", "C_freq", "G_freq", "A_freq", "Del_freq", "Ins_freq"], loc="right")
f.savefig(os.path.join(results_folder_path, f"{suptitle}.tiff"), dpi=300)

# save to disk
df_dna_forw_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.dna.forw.tsv"), sep="\t")
df_dna_rev_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.dna.rev.tsv"), sep="\t")
df_ribo_forw_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.ribo.forw.tsv"), sep="\t")
df_ribo_rev_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.ribo.rev.tsv"), sep="\t")

# ----- Start Computations ----------------#
# SITE 3 (rC)
site = 5004
start = site - 7
stop = site + 7
suptitle = f"DNA only vs Ribo1 site M13mp18:{site}"
print(f"[{datetime.now()}] Processing site {suptitle}.", flush=True)

# dna only forward
df_dna_forw_s1 = get_stats_from_region(dna_forw_path, "M13mp18", start, stop)
# dna only reverse
df_dna_rev_s1 = get_stats_from_region(dna_rev_path, "M13mp18", start, stop)
# Ribo forward
df_ribo_forw_s1 = get_stats_from_region(ribo_forw_path, "M13mp18", start, stop)
# Ribo reverse
df_ribo_rev_s1 = get_stats_from_region(ribo_rev_path, "M13mp18", start, stop)

# PLOTTING RESULTS
f, axes = plt.subplots(nrows=2, ncols=2, figsize=(15,7))
# plot results of dna only forward strand
df_dna_forw_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[0,0], legend=None,
                                color=["red", "blue", "orange", "green", "gray", "black"])
axes[0,0].set_title("DNA_only (+) strand")
axes[0,0].set_xticks([])
# plot results of ribo forward strand
df_ribo_forw_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[0,1], legend=None,
                                  color=["red", "blue", "orange", "green", "gray", "black"])
axes[0,1].set_title("Ribo (+) strand")
axes[0,1].set_xticks([])
# plot results for dna only reverse strand
df_dna_rev_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[1,0], legend=None,
                               color=["red", "blue", "orange", "green", "gray", "black"])
axes[1,0].set_title("DNA_only (-) strand")
# plot results for Ribo reverse strand
df_ribo_rev_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[1,1], legend=None, 
                                 color=["red", "blue", "orange", "green", "gray", "black"])
axes[1,1].set_title("Ribo (-) strand")
f.suptitle(suptitle)
axes[0,0].set_ylabel("Frequency")
axes[1,0].set_ylabel("Frequency")
axes[1,0].set_xlabel("Reference Position")
axes[1,1].set_xlabel("Reference Position")
f.legend(["T_freq", "C_freq", "G_freq", "A_freq", "Del_freq", "Ins_freq"], loc="right")
f.savefig(os.path.join(results_folder_path, f"{suptitle}.tiff"), dpi=300)

# save to disk
df_dna_forw_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.dna.forw.tsv"), sep="\t")
df_dna_rev_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.dna.rev.tsv"), sep="\t")
df_ribo_forw_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.ribo.forw.tsv"), sep="\t")
df_ribo_rev_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.ribo.rev.tsv"), sep="\t")

# ----- Start Computations ----------------#
# SITE 4 (U)
site = 5015
start = site - 7
stop = site + 7
suptitle = f"DNA only vs Ribo1 site M13mp18:{site}"
print(f"[{datetime.now()}] Processing site {suptitle}.", flush=True)

# dna only forward
df_dna_forw_s1 = get_stats_from_region(dna_forw_path, "M13mp18", start, stop)
# dna only reverse
df_dna_rev_s1 = get_stats_from_region(dna_rev_path, "M13mp18", start, stop)
# Ribo forward
df_ribo_forw_s1 = get_stats_from_region(ribo_forw_path, "M13mp18", start, stop)
# Ribo reverse
df_ribo_rev_s1 = get_stats_from_region(ribo_rev_path, "M13mp18", start, stop)

# PLOTTING RESULTS
f, axes = plt.subplots(nrows=2, ncols=2, figsize=(15,7))
# plot results of dna only forward strand
df_dna_forw_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[0,0], legend=None,
                                color=["red", "blue", "orange", "green", "gray", "black"])
axes[0,0].set_title("DNA_only (+) strand")
axes[0,0].set_xticks([])
# plot results of ribo forward strand
df_ribo_forw_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[0,1], legend=None,
                                  color=["red", "blue", "orange", "green", "gray", "black"])
axes[0,1].set_title("Ribo (+) strand")
axes[0,1].set_xticks([])
# plot results for dna only reverse strand
df_dna_rev_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[1,0], legend=None,
                               color=["red", "blue", "orange", "green", "gray", "black"])
axes[1,0].set_title("DNA_only (-) strand")
# plot results for Ribo reverse strand
df_ribo_rev_s1.iloc[:,-6:].plot(kind="bar", stacked=True, ax=axes[1,1], legend=None, 
                                 color=["red", "blue", "orange", "green", "gray", "black"])
axes[1,1].set_title("Ribo (-) strand")
f.suptitle(suptitle)
axes[0,0].set_ylabel("Frequency")
axes[1,0].set_ylabel("Frequency")
axes[1,0].set_xlabel("Reference Position")
axes[1,1].set_xlabel("Reference Position")
f.legend(["T_freq", "C_freq", "G_freq", "A_freq", "Del_freq", "Ins_freq"], loc="right")
f.savefig(os.path.join(results_folder_path, f"{suptitle}.tiff"), dpi=300)

# save to disk
df_dna_forw_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.dna.forw.tsv"), sep="\t")
df_dna_rev_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.dna.rev.tsv"), sep="\t")
df_ribo_forw_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.ribo.forw.tsv"), sep="\t")
df_ribo_rev_s1.to_csv(os.path.join(results_folder_path, f"{suptitle}.ribo.rev.tsv"), sep="\t")


###--- printing end of computations ---###
stop_time = datetime.now()
print(f"[{datetime.now()}] Computation Finished.", flush=True)
print(f"[{datetime.now()}] Elapsed time {stop_time-start_time}", flush=True)