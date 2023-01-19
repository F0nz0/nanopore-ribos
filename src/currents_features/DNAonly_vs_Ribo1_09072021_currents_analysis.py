# RIBO1 rep2 09/07/2021 analysis

# importing basic modules
import pandas as pd
import numpy as np
import os
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sn

# to plot better figures
plt.tight_layout()


# ----- defining basic paths -----#
dna_only_folder = "/lustre/bio_running/DNA_Ribo/run_dna_30062021"
ribo_folder = "/lustre/bio_running/DNA_Ribo/run_1xRIBO_09072021"
dna_only_eventalign_filename = "dna_only_wflag_wholeprimer.eventalign"
ribo_eventalign_filename = "run_1xRIBO_09072021_wflag_wholeprimer.eventalign"
dna_eventalign_path = os.path.join(dna_only_folder, "f5c_eval_wflags", dna_only_eventalign_filename)
ribo_eventalign_path = os.path.join(ribo_folder, "f5c_eval_wflags", ribo_eventalign_filename)
results_folder_path = "/lustre/bio_running/DNA_Ribo/results"


# ----- utils functions -----#
# def. aggregation function for samples
def agg_samples(x):
    '''
    Function to aggregate samples feature in the groupby operation
    merging all samples from the same event mapping onto a position
    into a unique list of float values 
    '''
    final = list( map(float, ",".join(x).split(",")))
    #final = ""
    #for e in x:
        #final += e
        #final += ","
        #final_list = final.rstrip(",").split(",")
        #final_map = list(map(float, final_list))
    return final


def is_rev_compl(seq1, seq2):
    '''
    Function to asses if two sequences are reverse complements.
    (for DNA)
    '''
    compl_dict = {"A":"T", "C":"G", "G":"C", "T":"A"}
    seq1_rev = seq1[::-1]
    seq1_rc = []
    for base in seq1_rev:
        seq1_rc.append(compl_dict[base])
    seq1_rc = "".join(seq1_rc)
    if seq1_rc == seq2:
        return True
    elif seq1_rc != seq2:
        return False
    

def collapse_eventalign(eventalign_filepath):
    '''
    Function to collapse on contig-position-ref_kmer-read_name the events from
    an eventalign file.
    '''
    df = pd.read_table(eventalign_filepath)
    # removing rows with NNNNNN model_kmer
    df = df[df["model_kmer"] != "NNNNNN"]
    df = df.groupby(["contig", "position", "reference_kmer", "read_name"]).agg({"samples":[agg_samples]})
    # resetting index and columns' names
    df = df.reset_index()
    df.columns = df.columns.droplevel(-1)
    df["event_level_mean"] = df["samples"].apply(np.mean)
    df["event_level_std"] = df["samples"].apply(np.std)
    df["dwell"] = df["samples"].apply(len).values
    df.drop("samples", axis=1, inplace=True)
    df.sort_values(["contig", "read_name"], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def eventalign_split_for_rev(eventalign_path):
    '''
    Custom function to split an eventalign file into forward and reverse tables modifying
    accordingly the strand field 
    + --> forward 
    +p --> palindromic reverse complement assigned as forward 
    - --> reverse
    '''
    # open input eventalign file
    eventalign = open(eventalign_path, "r")

    # open output file for forward and reverse eventalign events
    out_forward = open(eventalign_path + ".forw_strand", "wt")
    out_reverse = open(eventalign_path + ".rev_strand", "wt")

    rev_compl_count = 0
    rev_compl_notPal_count = 0
    for c, l in enumerate(eventalign):
        if l.startswith('contig'):
            out_forward.write(l)
            out_reverse.write(l)
            continue
        line = l.rstrip().split("\t")
        if is_rev_compl(line[2], line[9]):
            rev_compl_count += 1
            if line[2] != line[9]:
                rev_compl_notPal_count += 1
                line[4] = "-"
                out_reverse.write("\t".join(line)+"\n")
            else:
                line[4] = "+p"
                out_forward.write("\t".join(line)+"\n")
        else:
            line[4] = "+"
            out_forward.write("\t".join(line)+"\n")

    print(f"[{datetime.now()}] Total rows evaluated from the eventalign file:", c, flush=True)
    print(f"[{datetime.now()}] Total rows assigned as forward:", c - rev_compl_notPal_count, flush=True)
    print(f"[{datetime.now()}] Total reverse complement found:", rev_compl_count, flush=True)
    print(f"[{datetime.now()}] Total reverse complement found (palindromic excluded):", rev_compl_notPal_count, flush=True)

    eventalign.close()
    out_forward.close()
    out_reverse.close()

# ------------ START COMPUTATIONS --------------#
start_time = datetime.now()
print(f"[{datetime.now()}] Starting computations.", flush=True)

# splitting eventalign table for dna only run
print(f"[{datetime.now()}] Splitting eventalign file forw and rev from file {dna_eventalign_path}", flush=True)
eventalign_split_for_rev(dna_eventalign_path)

# splitting eventalign table for ribo run
print(f"[{datetime.now()}] Splitting eventalign file forw and rev from file {ribo_eventalign_path}", flush=True)
eventalign_split_for_rev(ribo_eventalign_path)

# create collapsed versions for each dataframes
# DNA ONLY
# forward strand
dna_forward_eventalign_filepath = dna_eventalign_path + ".forw_strand"
print(f"[{datetime.now()}] Collapsing events for {dna_forward_eventalign_filepath}", flush=True)
dna_forward_eventalign = collapse_eventalign(dna_forward_eventalign_filepath)

# reverse strand
dna_reverse_eventalign_filepath = dna_eventalign_path + ".rev_strand"
print(f"[{datetime.now()}] Collapsing events for {dna_reverse_eventalign_filepath}", flush=True)
dna_reverse_eventalign = collapse_eventalign(dna_reverse_eventalign_filepath)

# RIBO
# forward strand
ribo_forward_eventalign_filepath = ribo_eventalign_path + ".forw_strand"
print(f"[{datetime.now()}] Collapsing events for {ribo_forward_eventalign_filepath}", flush=True)
ribo_forward_eventalign = collapse_eventalign(ribo_forward_eventalign_filepath)

# reverse strand
ribo_reverse_eventalign_filepath = ribo_eventalign_path + ".rev_strand"
print(f"[{datetime.now()}] Collapsing events for {ribo_reverse_eventalign_filepath}", flush=True)
ribo_reverse_eventalign = collapse_eventalign(ribo_reverse_eventalign_filepath)


# PLOT RESULTS
# Violin Plots distributions
# FORWARD STRAND
# merging datasets
print(f"[{datetime.now()}] Concatenating forward currents.", flush=True)
df_forward_merged = pd.concat([dna_forward_eventalign, ribo_forward_eventalign], ignore_index=True)
df_forward_merged["position"] = df_forward_merged["position"] + 1 # change to 1-based coordinate system since eventaling is 0-based
df_forward_merged["y"] = ["DNAonly" for i in range(dna_forward_eventalign.shape[0])] + ["Ribo1" for i in range(ribo_forward_eventalign.shape[0])]

# make free some memory
del dna_forward_eventalign
del ribo_forward_eventalign

# REVERSE STRAND
# merging datasets
print(f"[{datetime.now()}] Concatenating reverse currents.", flush=True)
df_reverse_merged = pd.concat([dna_reverse_eventalign, ribo_reverse_eventalign], ignore_index=True)
df_reverse_merged["position"] = df_reverse_merged["position"] + 1 # change to 1-based coordinate system since eventaling is 0-based
df_reverse_merged["y"] = ["DNAonly" for i in range(dna_reverse_eventalign.shape[0])] + ["Ribo1" for i in range(ribo_reverse_eventalign.shape[0])]

# make free some memory
del dna_reverse_eventalign
del ribo_reverse_eventalign


# PLOTTING
print(f"[{datetime.now()}] Plotting results.", flush=True)

# SITE 1 (rA)
# defining sites and start and stop coordinates
site = 4985
start = site - 7
stop = site + 7

title = "DNAonly vs Ribo1 - rep2"
suptitle = f"M13mp18:{site}"

# defining subplots
f, axes = plt.subplots(nrows=2, ncols=1, figsize=(18, 9), sharex=True)

# plot figure for forward strand on axes 1
sn.violinplot(data=df_forward_merged.query(f"position >= {start} and position <= {stop}"), 
              x="position", 
              y="event_level_mean", 
              hue="y", 
              split=True,
              scale="width",
              ax=axes[0])
axes[0].set_title(title + ": Forward Strand")
axes[0].set_ylabel("Currents (pA)")
axes[0].legend(loc="upper right")

# plot figure for reverse strand on axes 2
sn.violinplot(data=df_reverse_merged.query(f"position >= {start} and position <= {stop}"), 
              x="position", 
              y="event_level_mean", 
              hue="y", 
              split=True,
              scale="width",
              ax=axes[1])
axes[1].set_title(title + ": Reverse Strand")
axes[1].set_ylabel("Currents (pA)")
axes[1].legend('', frameon=False)

f.suptitle(suptitle, size=15, x=0.515, y=0.945)
plt.savefig(os.path.join(results_folder_path, f"{title}-{suptitle}.tiff"), dpi=300)
print(f"[{datetime.now()}] Figure Saved.", flush=True)

# SITE 2 (rG)
# defining sites and start and stop coordinates
site = 4997
start = site - 7
stop = site + 7

title = "DNAonly vs Ribo1 - rep2"
suptitle = f"M13mp18:{site}"

# defining subplots
f, axes = plt.subplots(nrows=2, ncols=1, figsize=(18, 9), sharex=True)

# plot figure for forward strand on axes 1
sn.violinplot(data=df_forward_merged.query(f"position >= {start} and position <= {stop}"), 
              x="position", 
              y="event_level_mean", 
              hue="y", 
              split=True,
              scale="width",
              ax=axes[0])
axes[0].set_title(title + ": Forward Strand")
axes[0].set_ylabel("Currents (pA)")
axes[0].legend(loc="upper right")

# plot figure for reverse strand on axes 2
sn.violinplot(data=df_reverse_merged.query(f"position >= {start} and position <= {stop}"), 
              x="position", 
              y="event_level_mean", 
              hue="y", 
              split=True,
              scale="width",
              ax=axes[1])
axes[1].set_title(title + ": Reverse Strand")
axes[1].set_ylabel("Currents (pA)")
axes[1].legend('', frameon=False)

f.suptitle(suptitle, size=15, x=0.515, y=0.945)
plt.savefig(os.path.join(results_folder_path, f"{title}-{suptitle}.tiff"), dpi=300)
print(f"[{datetime.now()}] Figure Saved.", flush=True)

# SITE 3 (rC)
# defining sites and start and stop coordinates
site = 5004
start = site - 7
stop = site + 7

title = "DNAonly vs Ribo1 - rep2"
suptitle = f"M13mp18:{site}"

# defining subplots
f, axes = plt.subplots(nrows=2, ncols=1, figsize=(18, 9), sharex=True)

# plot figure for forward strand on axes 1
sn.violinplot(data=df_forward_merged.query(f"position >= {start} and position <= {stop}"), 
              x="position", 
              y="event_level_mean", 
              hue="y", 
              split=True,
              scale="width",
              ax=axes[0])
axes[0].set_title(title + ": Forward Strand")
axes[0].set_ylabel("Currents (pA)")
axes[0].legend(loc="upper right")

# plot figure for reverse strand on axes 2
sn.violinplot(data=df_reverse_merged.query(f"position >= {start} and position <= {stop}"), 
              x="position", 
              y="event_level_mean", 
              hue="y", 
              split=True,
              scale="width",
              ax=axes[1])
axes[1].set_title(title + ": Reverse Strand")
axes[1].set_ylabel("Currents (pA)")
axes[1].legend('', frameon=False)

f.suptitle(suptitle, size=15, x=0.515, y=0.945)
plt.savefig(os.path.join(results_folder_path, f"{title}-{suptitle}.tiff"), dpi=300)
print(f"[{datetime.now()}] Figure Saved.", flush=True)

# SITE 3 (U)
# defining sites and start and stop coordinates
site = 5015
start = site - 7
stop = site + 7

title = "DNAonly vs Ribo1 - rep2"
suptitle = f"M13mp18:{site}"

# defining subplots
f, axes = plt.subplots(nrows=2, ncols=1, figsize=(18, 9), sharex=True)

# plot figure for forward strand on axes 1
sn.violinplot(data=df_forward_merged.query(f"position >= {start} and position <= {stop}"), 
              x="position", 
              y="event_level_mean", 
              hue="y", 
              split=True,
              scale="width",
              ax=axes[0])
axes[0].set_title(title + ": Forward Strand")
axes[0].set_ylabel("Currents (pA)")
axes[0].legend(loc="upper right")

# plot figure for reverse strand on axes 2
sn.violinplot(data=df_reverse_merged.query(f"position >= {start} and position <= {stop}"), 
              x="position", 
              y="event_level_mean", 
              hue="y", 
              split=True,
              scale="width",
              ax=axes[1])
axes[1].set_title(title + ": Reverse Strand")
axes[1].set_ylabel("Currents (pA)")
axes[1].legend('', frameon=False)

f.suptitle(suptitle, size=15, x=0.515, y=0.945)
plt.savefig(os.path.join(results_folder_path, f"{title}-{suptitle}.tiff"), dpi=300)
print(f"[{datetime.now()}] Figure Saved.", flush=True)

# printing end of computations
stop_time = datetime.now()
print(f"[{datetime.now()}] Computation Finished.", flush=True)
print(f"[{datetime.now()}] Elapsed time {stop_time-start_time}", flush=True)