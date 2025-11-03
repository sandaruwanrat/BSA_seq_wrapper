#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import re
from datetime import datetime
from pathlib import Path
import time
import pandas as pd
import seaborn as sns
from statsmodels.stats.proportion import proportion_confint
import matplotlib
import matplotlib.pyplot as plt


# ===========================================================

datestring = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")


# ===================== Argument parser ======================================

# input argument parser setup
parser = argparse.ArgumentParser(
    description="Plot Amplicon sequencing data",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

# arguments for default mode
parser.add_argument('-R', '--files', default="", help="Directory containing Ref fasta")
parser.add_argument('-D', '--matpath', default="", help="Directory of output matrix")

# optional arguments
parser.add_argument('-v', '--verbose', action='store_true',
                    help="Verbose about parameters and steps (subset of log file).")

# Parse arguments
args = parser.parse_args()


# ===================== Assign variables ======================================

dir_ref = args.files
dir_out = args.matpath
verbose = args.verbose


#_______________________________________________________________________
# ===================== Functions ======================================
#_______________________________________________________________________



# ===================== get unique amplicons from ref fasta ======================================

def get_amplicon_list(ref_dir: Path):
    amplicon_list = []
    
    # iterate over all fasta files in directory
    #later add a funtion to look for multiple file type .fa or fasta
 
    for fa_file in Path(ref_dir).glob("*.fa"):
        #print(fa_file)
        print("fa file: {}".format(fa_file))
        with open(fa_file, 'r') as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    match = re.match(r'^>(\S+)', line)
                    print(match)
                    if match:
                        amplicon_list.append(match.group(1))

    # remove duplicates while preserving order
    seen = set()
    unq_amplicons = [x for x in amplicon_list if not (x in seen or seen.add(x))]
    print("unq_amplicons : {}".format(unq_amplicons))
    return unq_amplicons


# ===================== Cal wilson CI ======================================

def cal_wilson_CI(dir_out_pth: Path):
    
    out_mat_file_pth=os.path.join(dir_out_pth, "output_matrix.txt")
    
    mat_df=pd.read_csv(out_mat_file_pth,header=0,sep='\t')
   # mat_df = pd.read_csv(out_mat_file_pth, sep='\t', names=["amplicon","context","percentage","reads","strand"], header=None)

    
    mat_df.rename(columns = {"reads": "totalC"}, 
            inplace = True)
    mat_df.loc[:,'Per']=mat_df["percentage"]*100
    
    df_filtered = mat_df[mat_df["strand"] == "+"].copy()
    df_filtered.loc[:,"count"] = (df_filtered["Per"] / 100) * df_filtered["totalC"]
    
    # Compute Wilson 95% CI
    lower, upper = proportion_confint(
        count=df_filtered['count'],
        nobs=df_filtered['totalC'],
        alpha=0.05,
        method='wilson'
    )

    # Add to DataFrame
    df_filtered['lower'] = lower * 100
    df_filtered['upper'] = upper * 100

    # Compute offsets (for error bars)
    df_filtered['plus'] = df_filtered['upper'] - df_filtered['Per']
    df_filtered['minus'] = df_filtered['Per'] - df_filtered['lower']
    
    
    return df_filtered


# ===================== plot amplicon specific methylation percentage======================================
def plot_amplicon(ref_dir: Path,out_mat_pth:Path):
    """
    Generate bar plots for each unique amplicon and save as PDF.

    Parameters
    ----------
    ref_dir : Path
        Directory containing reference fasta
    plot_out_path : Path
        Directory where PDF plots will be saved
    
    """
    #cal CI for output
    plot_df=cal_wilson_CI(out_mat_pth)
    plot_df.loc[:,["Amp"]] =  plot_df["amplicon"].str.rsplit("_", n=1, expand=True)[0]
    plot_df.loc[:,["Smp"]] =  plot_df["amplicon"].str.rsplit("_", n=1, expand=True)[1]
    
    #get the unique amplicons from fasta headers
    unique_amps = get_amplicon_list(ref_dir)
    print("unique_amps: {}".format(unique_amps))
    
    
    #plotting
    #color palette
    
    plot_out_path = Path(out_mat_pth)
    #plot_out_path.mkdir(parents=True, exist_ok=True)
    
    palette = ["#EE3624", "#42A9DF", "#00A850"]

    for amp in unique_amps:
        # Subset data by amplicon
        subset =  plot_df[plot_df["Amp"] == amp]

        plt.figure(figsize=(10, 6))
        
        #python >3.7 errorbar=None
        sns.barplot(
            data=subset,
            x="amplicon",
            y="Per",
            hue="context",
            palette=palette,
            ci=None
        )

        # Add error bars manually (Per ± plus/minus)
        for i, (name, group) in enumerate(subset.groupby("context")):
            # Calculate x positions for grouped bars
            x_positions = range(len(group))
            offsets = [-0.27, 0, 0.27]  # adjust spacing manually for clarity if needed
            plt.errorbar(
                x=[x + offsets[i] for x in range(len(group))],
                y=group["Per"],
                yerr=[group["minus"], group["plus"]],
                fmt="none",
                capsize=4,
                color="black",
                elinewidth=1,
            )

        # Add horizontal line (y = 0)
        plt.axhline(0, color="black", linewidth=1)
            
        plt.ylim(0, 100)
        plt.ylabel("Percentage")
        plt.xlabel(f"Amplicon : {amp}")
        plt.axhline(0, color="black", linewidth=1)
        plt.xticks(rotation=90, ha="center")
        plt.legend(title="Context")

        plt.tight_layout()

        # Save to PDF (one file per amplicon)
        pdf_path = plot_out_path / f"{amp}_plot.pdf"
        plt.savefig(pdf_path, format="pdf", bbox_inches="tight")

        plt.close()  

    return True



#___________________________________________________________________
# ===================== main ======================================
#___________________________________________________________________

def main():
    print("Output directory: {}".format(dir_out))
    print("Ref directory: {}".format(dir_ref))
    print("\n")
    

    # directory paths
    dir_out_pth = Path(dir_out).resolve()
    print("output path: {}".format(dir_out_pth))
    dir_ref_pth = Path(dir_ref).resolve()
    print("Ref path: {}".format(dir_ref_pth))
  
    
    # Chek output directory
    if dir_out_pth.exists():
        if args.verbose:
            print(f"[INFO] outmat directory exists: {dir_out_pth}")
    else:
        if args.verbose:
            print(f"[INFO] wrong outmat directory: {dir_out_pth}")
        
        
    # Initialize logfiles
    datestring_log = datetime.now().strftime("%I:%M:%S")
    
    log_file = Path(dir_out_pth) / "plot_amp.log"
    
    with open(log_file, "w") as plt_log:
        if verbose:
            print(f"[INFO] Starting run at {datestring_log}", file=sys.stderr)
        print(f"Starting run at {datestring_log }", file=plt_log)
        
        log_msg = (
            f"Reference_dir = {dir_ref_pth}\n"
            f"Output_dir = {dir_out_pth}\n"
        )
        print(log_msg, file=plt_log)
        if verbose:
            print(log_msg, file=sys.stderr)

        
        #plotting
        print("\tPlotting output matix results", file=plt_log)
        if verbose:
            print("\tPlotting output matix results", file=sys.stderr)
        
        plot_amplicon(ref_dir=dir_ref_pth,out_mat_pth=dir_out_pth)
        
        #write_log(fh_log, "\t" + " ".join(cmd_BI), verbose)
        
    
        
          
     


if __name__ == "__main__":
    
    start_time = time.time()
    # Missing argument checks
    if not dir_ref:
        print("Error >>>>>>>>>>>> Directory containing reference files not provided <<<<<<<<<<<<<<<<<")
        
    elif not dir_out:
        print("Error >>>>>>>>>>>> output mat dir not found <<<<<<<<<<<<<<<<<")
        
    

    
    # Path existence checks
    if not os.path.isdir(dir_ref):
        print("Error >>>>>>>>>>>> Reference directory does not exist: {} <<<<<<<<<<<<<<<<<".format(dir_ref))
        sys.exit(1)
    if not os.path.isdir(dir_out):
        print("Error >>>>>>>>>>>> output mat dir does not exist: {} <<<<<<<<<<<<<<<<<".format(dir_out))
        sys.exit(1)
    
    
    main()
    print("Amplicon Plot completed")

    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"Elapsed time: {elapsed_time:.4f} seconds")




