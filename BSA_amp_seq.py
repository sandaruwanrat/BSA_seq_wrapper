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
import matplotlib as plt



# Define timestamp for output folder
datestring = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")



# input argument parser setup
parser = argparse.ArgumentParser(
    description="Analyze Bisulfite Amplicon sequencing reads - map reads using methylpy, generate allC files, count methylation percentage, coverage,stats and plots data",
    formatter_class=argparse.RawDescriptionHelpFormatter
)

# arguments for default mode
parser.add_argument('-S', '--sequences', default="", help="Directory containing the sequencing reads files (e.g., /cluster/pixstor/slotkinr-lab/sandaruwan/BSA_round18/KS-SAE-Amplicon-R18).")
parser.add_argument('-F', '--files', default="", help="Directory containing the Reference files (e.g., ~/BSA_round18/REF_Files/).")
parser.add_argument('-M', '--master', default="", help="Master file containing the sample names and the reference fasta file names.")

# optional arguments

parser.add_argument('-o', '--output', default=f"output_{datestring}",
                    help="Name of the directory to store output files (default: 'output_CURRENTTIME')")
parser.add_argument('-n', '--numberOfCHSites', type=int, default=5,
                    help="Minimum number of CH sites for mCH level filter (default: 5).")
parser.add_argument('-t', '--threshold', type=float, default=0.95,
                    help="Threshold of mCH percentage to filter reads (e.g., 96.78 filters >96.78%% CH sites). Use 101 to disable (default: 95).")
parser.add_argument('-p', '--processors', type=int, default=2,
                    help="Number of processors for methylpy (default: 2).")
parser.add_argument('-v', '--verbose', action='store_true',
                    help="Verbose about parameters and steps (subset of log file).")

# Parse arguments
args = parser.parse_args()

# Assign variables
dir_reads = args.sequences
dir_ref = args.files
dir_out = args.output
master_file = args.master
filT = args.threshold
num_proc = args.processors
siteCH = args.numberOfCHSites
verbose = args.verbose



# Help function
def print_help():
    print(f"""
#### {os.path.basename(sys.argv[0])} ####

#note:
##later plotting is used, need to add a input asking for output plot using ifelse

# This script is used to analyze Bisulfite Amplicon sequencing reads - map reads using methylpy, generate allC files, count methylation percentage, coverage and stats.

Usage: {os.path.basename(sys.argv[0])} -S <dir> -F <dir> -M <file> OPTIONS

Required Inputs (default mode OR Re-Run mode):

    For default mode:
        -S|--sequences <dir>    Directory containing the sequencing reads files (e.g., /cluster/yourname/BSA_round18/KS-SAE-Amplicon-R1).
        -F|--files <dir>        Directory containing the reference files (e.g., ~/BSA_round18/REF_Files/).
        -M|--master <file>      Master file containing the sample names and the reference fasta file names.

    Optional Inputs:
        -o|--output <dir>       Name of the directory to create and store the output files in (e.g., output1). (default: 'output_CURRENTTIME')
        -n|--numberOfCHSites <int>  Minimum number of CH sites for mCH level filter to operate. (default: 5)
        -t|--threshold <float>  Threshold of mCH percentage above which reads will be filtered out. (default: 95)
        -p|--processors <int>   Number of processors to be used for methylpy. (default: 2)
        -v|--verbose            Verbose about parameters and steps.
        --help                  Print this help
""")
    sys.exit(1)


def timestamp_fn():
    return datetime.now().strftime("%Y_%m_%d_%H_%M_%S")


def ensure_dir_fn(path: Path, verbose: bool = False):
    if path.exists():
        if verbose:
            print(f"[INFO] Output directory already exists: {path}")
    else:
        if verbose:
            print(f"[INFO] Creating output directory: {path}")
        path.mkdir(parents=True, exist_ok=True)


def open_logfile_fn(outdir: Path):
    logfile = outdir / "amplicon_seq.log"
    return logfile.open("w")


def write_log(fh, msg: str, verbose=True):
    print(msg, file=fh)
    if verbose:
        print(msg, file=sys.stderr)
    
# --- Setup Functions ---
def initialize_output_files(outdir: Path):
    """Prepare empty results tables"""
    out_mat_file = outdir / "output_matrix.txt"
    with open(out_mat_file, "w") as fh:
        print("amplicon\tcontext\tpercentage\tsites\treads\tstrand", file=fh)

    stat_file = outdir / "stats_conversion.txt"
    with open(stat_file, "w") as fh:
        print("genotype\tconversion rate", file=fh)

    cov_file = outdir / "stats_coverage.txt"
    with open(cov_file, "w") as fh:
        print("amplicon\tstrand\t#totalCytosines\t#reportedCytosines\t%reportedCytosines\ttotalReads\tcoverage", file=fh)


def parse_master(master_file: str):
    """Yield (sample, fasta) from master file"""
    with open(master_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            match = re.match(r'(\S+)\t(\S+)\.fa', line)
            if match:
                yield match.group(1), match.group(2)


# def remove_bam(fh_log,bamfile_pth:Path):
#     """Remove unfiltered bamfiles"""
#     rm_bm_cmd = f"rm {bamfile_pth}"
#     print(f"\t{rm_bm_cmd}", file=fh_log)
#     subprocess.run(rm_bm_cmd, shell=True)
    
#     bam_index = f"{bamfile_pth}.bai"
#     if os.path.exists(bam_index):
#         rm_ind_cmd = f"rm {bam_index}"
#         print(f"\t{rm_ind_cmd}", file=fh_log)
#         subprocess.run(rm_ind_cmd, shell=True)
        
#     return True

def remove_bam(fh_log, bamfile_pth: Path):
    """Remove BAM file and its index if they exist."""
    if not bamfile_pth.exists():
        log_msg = f"\tNo BAM file found at {bamfile_pth}, skipping removal."
        print(log_msg, file=fh_log)
        return False

    # Remove BAM file
    rm_bm_cmd = f"rm {bamfile_pth}"
    print(f"\t{rm_bm_cmd}", file=fh_log)
    subprocess.run(rm_bm_cmd, shell=True)

    # Remove BAM index 
    bam_index = bamfile_pth.with_suffix(bamfile_pth.suffix + ".bai")
    if bam_index.exists():
        rm_ind_cmd = f"rm {bam_index}"
        print(f"\t{rm_ind_cmd}", file=fh_log)
        subprocess.run(rm_ind_cmd, shell=True)

    return True


# methylpy Commands
def run_build_reference(name_fasta, fasta_file, dir_ref_pth, num_proc, fh_log, verbose=True):
    cmd_BI = [
        "methylpy", "build-reference",
        "--num-procs", str(num_proc),
        "--input-files", str(fasta_file),
        "--output-prefix", f"{dir_ref_pth}/{name_fasta}_methylpy",
        "--aligner", "bowtie2"
    ]
    write_log(fh_log, "\t" + " ".join(cmd_BI), verbose)
    result = subprocess.run(cmd_BI)
    print("Completed indexing: {}".format(name_fasta))
    
    
    if result.returncode != 0:
        err_msg = f"Error building reference index for {name_fasta}"
        print(err_msg, file=sys.stderr)
        print(err_msg, file=fh_log)
        raise RuntimeError(err_msg)
    
    return True

    

#data processing for methylation calling, using single end pipeline
def run_single_end_pipeline(name_sample,bam_path, fastq_path, fwd_ref_file, rv_ref_file,
                            ref_fasta_pth, num_proc, dir_out_pth, mpy_er_path,
                            fh_log, verbose=True):
    # Build command as a string
    cmd_MP = (
        f"methylpy single-end-pipeline "
        f"--read-files {fastq_path} "
        f"--sample {name_sample} "
        f"--forward-ref {fwd_ref_file} "
        f"--reverse-ref {rv_ref_file} "
        f"--ref-fasta {ref_fasta_pth} "
        f"--num-procs {num_proc} "
        f"--min-cov 3 "
        f"--remove-chr-prefix False "
        f"--pbat False "
        f"--adapter-seq CTGTCTCTTATACACATCT "
        f"--compress-output False "
        f"--path-to-output {dir_out_pth}"
    )

    # Log command
    write_log(fh_log, "\t" + cmd_MP, verbose)

    # Run with shell=True since command is a string
    with open(mpy_er_path, "w") as er_file:
        result = subprocess.run(cmd_MP, shell=True, stdout=er_file, stderr=subprocess.STDOUT)

    if result.returncode != 0:
        err_msg = f"Error running methylpy for {name_sample}\nCheck {mpy_er_path}"
        write_log(fh_log, err_msg, verbose)
        print(err_msg, file=sys.stderr)
        raise RuntimeError(err_msg)
    # else:
    #     log_msg = f"\tFinished running methylpy for {name_sample}"
    #     write_log(fh_log, log_msg, verbose)
    #     if verbose:
    #         print(log_msg, file=sys.stderr)
            
    
    if not bam_path.exists():
        err_msg = f"BAM file not found for {name_sample} at {bam_path}"
        write_log(fh_log, err_msg, verbose)
        print(err_msg, file=sys.stderr)
        raise FileNotFoundError(err_msg)

    else:
        log_msg = f"\tFinished running methylpy for {name_sample} "
        write_log(fh_log, log_msg, verbose)
        if verbose:
            print(log_msg, file=sys.stderr)

    return True

def bqf_fn(name_sample,bam_path,filtered_bam_path, fasta_file,siteCH ,filT, filtering_er_path,
                            fh_log, verbose=True):
    
    
    cmd_BQ = (f"methylpy bam-quality-filter "
            f"--input-file {bam_path} "
            f"--output-file {filtered_bam_path} "
            f"--ref-fasta {fasta_file} "
            f"--min-num-ch {siteCH} --max-mch-level {filT} "
        )
    
    print(f"\t{cmd_BQ}", file=fh_log)
    
    #result = subprocess.run(cmd, shell=True)
    with open(filtering_er_path, "w") as filt_er_file:
            result_BQ = subprocess.run(cmd_BQ, shell=True, stdout=filt_er_file, stderr=subprocess.STDOUT)
    
    
    if result_BQ.returncode != 0:
        err_msg = f"Error filtering reads for {name_sample} \n check the filter_er file"
        print(err_msg, file=sys.stderr)
        print(err_msg, file=fh_log)
        #raise RuntimeError(f"Error filtering reads for {name_sample}, see {filtering_er_path}")
                
    else:
        log_msg = f"\tFinished filtering reads for {name_sample}"
        print(log_msg, file=fh_log)
        if verbose:
            print(log_msg, file=sys.stderr)

    
    return True

#methylation calling function

def call_methylation_state(name_sample,filtered_bam,
                            ref_fasta_pth, num_proc, dir_out_pth, cms_er_path,
                            fh_log, verbose=True):
    """Function to run methylation calling"""
    
    if not filtered_bam.exists():
        log_msg = f"\tNo filtered BAM file found at {filtered_bam}"
        print(log_msg, file=fh_log)
        raise FileNotFoundError(log_msg)
    
    cmd_CMS = (f"methylpy call-methylation-state "
        f"--input-file {filtered_bam} "
        f"--sample {name_sample} "
        f"--ref-fasta {ref_fasta_pth} "
        f"--num-procs {num_proc} "
        f"--min-cov 3 "
        f"--remove-chr-prefix FALSE "
        f"--compress-output FALSE --paired-end FALSE "
        f"--path-to-output {dir_out_pth} "
        )
    
    print(f"\t{cmd_CMS}", file=fh_log)
    
    
    with open(cms_er_path, "w") as cms_er_file:
            result_CMS = subprocess.run(cmd_CMS, shell=True, stdout=cms_er_file, stderr=subprocess.STDOUT)
    
    
    if result_CMS.returncode != 0:
        err_msg = f"call-methylation-state for {name_sample} \n check the cms_er file"
        print(err_msg, file=sys.stderr)
        print(err_msg, file=fh_log)
        raise RuntimeError(f"Error filtering reads for {name_sample}, check {cms_er_path}")
        #continue        
    else:
        log_msg = f"\tcall-methylation-state for {name_sample}"
        print(log_msg, file=fh_log)
        if verbose:
            print(log_msg, file=sys.stderr)
    
    return True

def check_reference(fasta_file, fh_log):
    """Check reference fasta file exists."""
    print(f"\tReference fasta file is {fasta_file}", file=fh_log)
    if not os.path.exists(fasta_file):
        print(f"Error: Reference fasta file not found: {fasta_file}", file=sys.stderr)
        return False
    return True

def get_bed_name(fasta_file):
    """Extract amplicon name(s) from fasta headers."""
    bed_names = []
    with open(fasta_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                match = re.match(r'^>(\S+)', line)
                if match:
                    bed_names.append(match.group(1))
    return bed_names

def cal_cytosine_cov(bed_file,ref_fasta,fh_log,verbose=True):
    """Count total C and G from ref fasta files for each amplicon region"""
    
    # plus strand Cs
    cmd_c = f"bedtools getfasta -fi {ref_fasta} -bed {bed_file} -fo - | grep -v '>' | tr -cd C | wc -c"
    print(f"\t{cmd_c}", file=fh_log)
    #for python 3.7 and above
    #result = subprocess.run(cmd_c, shell=True, capture_output=True, text=True)
    result = subprocess.run(cmd_c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    plus_total_c = int(result.stdout.strip()) if result.stdout.strip() else 0

    # minus strand Gs
    cmd_g = f"bedtools getfasta -fi {ref_fasta} -bed {bed_file} -fo - | grep -v '>' | tr -cd G | wc -c"
    print(f"\t{cmd_g}", file=fh_log)
    #for python 3.7 and above
    #result = subprocess.run(cmd_g, shell=True, capture_output=True, text=True)
    result = subprocess.run(cmd_g, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    minus_total_c = int(result.stdout.strip()) if result.stdout.strip() else 0
    
    return plus_total_c, minus_total_c


def process_strand(allc_file, bed_file, strand, name_bed, name_sample, fh_log):
    """Process methylation contexts for one strand and return stats list."""
    
    contexts = ['CG[ATCG]', 'C[ATC]G', 'C[ATC][ATC]']
    context_names = ['CG', 'CHG', 'CHH']
    outputs = []

    for context, context_name in zip(contexts, context_names):
        cmd = (f"awk '$3==\"{strand}\"{{print $1\"\\t\"$2-1\"\\t\"$2\"\\t\"NR\"\\t.\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6}}' "
               f"{allc_file} | "
               f"bedtools intersect -wa -a - -b {bed_file} | "
               f"awk '$7~/{context}/{{sum1=sum1+$8;sum2=sum2+$9;sum3=sum3+1}}"
               f"END{{if(sum2>0) print \"{name_bed}_{name_sample}\\t{context_name}\\t\"sum1/sum2\"\\t\"sum3\"\\t\"sum2\"\\t{strand}\"; "
               f"else print \"{name_bed}_{name_sample}\\t{context_name}\\t0\\t0\\t0\\t{strand}\"}}' 2>/dev/null")

        print(f"\t{cmd}", file=fh_log)
        #for python 3.7
        #result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        out = result.stdout.strip()

        if not out:
            out = f"{name_bed}_{name_sample}\t{context_name}\t0\t0\t0\t{strand}"

        outputs.append(out)

    return outputs       
        
        

def cal_cov_stats(fasta_file,name_sample,dir_ref_pth,dir_out_pth,fh_log,verbose=True):
    
    """runner for one sample and reference."""
    #fasta_file = f"{dir_ref}{name_fasta}.fa"
    if not check_reference(fasta_file, fh_log):
        return
    
    bed_names = get_bed_name(fasta_file)
    print("bed regions: {}".format(bed_names))
    allc_file = Path(dir_out_pth)/f"allc_{name_sample}.tsv"
    
    
    if not os.path.exists(allc_file):
        print(f"Warning: allc file not found: {allc_file}", file=sys.stderr)
        raise RuntimeError(f"allc file not found {allc_file}")
    
    for name_bed in bed_names:
        bed_file = Path(dir_ref_pth)/f"{name_bed}.bed"
        #process_amplicon(name_bed, name_sample, fasta_file, bed_file, allc_file, dir_out, fh_log)
        if not os.path.exists(bed_file):
            print(f"Warning: BED file not found: {bed_file}", file=sys.stderr)
            raise RuntimeError(f"bedfile file {bed_file} not found for {name_sample}")
        
        plus_total_c, minus_total_c = cal_cytosine_cov(ref_fasta=fasta_file, bed_file=bed_file,fh_log= fh_log)
        
        # process both strands
        plus_outs = process_strand(allc_file, bed_file, "+", name_bed, name_sample, fh_log)
        minus_outs = process_strand(allc_file, bed_file, "-", name_bed, name_sample, fh_log)
        
        # write output_matrix
        with open(os.path.join(dir_out, "output_matrix.txt"), "a") as fh_matrix:
            for line in plus_outs + minus_outs:
                print(line, file=fh_matrix)
        
        
        # coverage stats (simple sums)
        plus_reported_c = sum(int(l.split("\t")[3]) for l in plus_outs if l.split("\t")[3].isdigit())
        plus_total_r = sum(int(l.split("\t")[4]) for l in plus_outs if l.split("\t")[4].isdigit())
        minus_reported_c = sum(int(l.split("\t")[3]) for l in minus_outs if l.split("\t")[3].isdigit())
        minus_total_r = sum(int(l.split("\t")[4]) for l in minus_outs if l.split("\t")[4].isdigit())

        plus_per_c = f"{plus_reported_c * 100 / plus_total_c:.2f}" if plus_total_c > 0 else "0.00"
        plus_cov = f"{plus_total_r / plus_reported_c:.2f}" if plus_reported_c > 0 else "0.00"
        minus_per_c = f"{minus_reported_c * 100 / minus_total_c:.2f}" if minus_total_c > 0 else "0.00"
        minus_cov = f"{minus_total_r / minus_reported_c:.2f}" if minus_reported_c > 0 else "0.00"

        amplicon_name = f"{name_bed}_{name_sample}"
        with open(os.path.join(dir_out, "stats_coverage.txt"), "a") as fh_stats:
            print(f"{amplicon_name}\tplus\t{plus_total_c}\t{plus_reported_c}\t{plus_per_c}\t{plus_total_r}\t{plus_cov}", file=fh_stats)
            print(f"{amplicon_name}\tminus\t{minus_total_c}\t{minus_reported_c}\t{minus_per_c}\t{minus_total_r}\t{minus_cov}", file=fh_stats)

            
    
    return True




def get_amplicon_list(ref_dir: Path):
    amplicon_list = []

    # iterate over all fasta files in directory
    for fa_file in Path(ref_dir).glob("*.fasta"):
        with open(fa_file, 'r') as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    match = re.match(r'^>(\S+)', line)
                    if match:
                        amplicon_list.append(match.group(1))

    # remove duplicates while preserving order
    seen = set()
    unq_amplicons = [x for x in amplicon_list if not (x in seen or seen.add(x))]

    return unq_amplicons


def cal_wilson_CI(dir_out_pth: Path):
    
    out_mat_file_pth=os.path.join(dir_out_pth, "output_matrix.txt")
    
    mat_df=pd.read_csv(out_mat_file_pth,header=0,sep='\t')
    
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

def plot_amplicon(ref_dir: Path,dir_out_pth:Path):
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
    plot_df=cal_wilson_CI(dir_out_pth)
    plot_df.loc[:,["Amp"]] =  plot_df["amplicon"].str.rsplit("_", n=1, expand=True)[0]
    plot_df.loc[:,["Smp"]] =  plot_df["amplicon"].str.rsplit("_", n=1, expand=True)[1]
    
    #get the unique amplicons from fasta headers
    unique_amps = get_amplicon_list(ref_dir)
    palette = ["#EE3624", "#42A9DF", "#00A850"]

    plot_out_path = Path(plot_out_path)
    plot_out_path.mkdir(parents=True, exist_ok=True)

    for amp in unique_amps:
        # Subset data by amplicon
        subset =  plot_df[plot_df["Amp"] == amp]

        plt.figure(figsize=(10, 6))
        sns.barplot(
            data=subset,
            x="amplicon",
            y="Per",
            hue="context",
            palette=palette,
            errorbar=None
        )

        # Add error bars manually (Per ± plus/minus)
        for i, (name, group) in enumerate(test_df.groupby("context")):
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
        pdf_path = dir_out_pth / f"{amp}_plot.pdf"
        plt.savefig(pdf_path, format="pdf", bbox_inches="tight")

        plt.close()  

    return True



#### main ######

def main():
    print("starting .....")
    print("date_time: {}".format(datestring))
    print("Reads directory: {}".format(dir_reads))
    print("Reference directory: {}".format(dir_ref))
    print("Output directory: {}".format(dir_out))
    print("Master file: {}".format(master_file))
    print("Filter threshold: {}".format(filT))
    print("Minimum number of CH sites: {}".format(siteCH))
    print("\n")
    
    
    #from perl may not need to be defined for now
    #my ($stats_file,$cmd,$index_file,$bam_file);
    
    
    # Initialize variables
    plusTotalC, plusReportedC, plusTotalR = 0, 0, 0
    minusTotalC, minusReportedC, minusTotalR = 0, 0, 0
    
    stats_file = None
    cmd = None
    index_file = None
    bam_file = None

    # Normalize directory paths
    dir_out_pth = Path(dir_out).resolve()
    print("output path: {}".format(dir_out_pth))
    dir_ref_pth = Path(dir_ref).resolve()
    print("Ref path: {}".format(dir_ref_pth))
    dir_reads_pth = Path(dir_reads).resolve()
    print("Input reads path: {}".format(dir_reads_pth))
    
    filT_converted = filT / 100  # Convert methylation percentage to decimal
    
    # Chek/create output directory
    if dir_out_pth.exists():
        if args.verbose:
            print(f"[INFO] Output directory already exists: {dir_out_pth}")
    else:
        if args.verbose:
            print(f"[INFO] Creating output directory: {dir_out_pth}")
        dir_out_pth.mkdir(parents=True, exist_ok=True)
        
        
    # Initialize logfiles
    datestring_log = datetime.now().strftime("%I:%M:%S")
    
    log_file = Path(dir_out_pth) / "amplicon_seq.log"
    
    with open(log_file, "w") as fh_log:
        if verbose:
            print(f"[INFO] Starting run at {datestring_log}", file=sys.stderr)
        print(f"Starting run at {datestring_log }", file=fh_log)
        
        log_msg = (
            f"Input_reads_dir = {dir_reads_pth}\n"
            f"Reference_dir = {dir_ref_pth}\n"
            f"Output_dir = {dir_out_pth}\n"
            f"Master file = {master_file}\n"
            f"Minimum number of CH sites for filtering = {siteCH}\n"
            f"Threshold for filtering reads = {filT}\n"
            f"Number of processors = {num_proc}"
        )
        print(log_msg, file=fh_log)
        if verbose:
            print(log_msg, file=sys.stderr)


        # Initialize output files
        initialize_output_files(dir_out_pth)
        
        

        # Process master file
        try:
            with open(master_file, 'r') as f_mas:
                for line in f_mas:
                    line = line.strip()
                    print("master file line: {}".format(line))
                    if not line:
                        continue
                        
                    match = re.match(r'(\S+)\t(\S+)\.fa', line)
                    if not match:
                        continue
                        
                    name_sample = match.group(1)
                    print("sample name: {}".format(name_sample))
                    name_fasta = match.group(2)
                    print("fasta_name: {}".format(name_fasta))
                    
                    print(f"Processing {name_sample}", file=fh_log)
                    if verbose:
                        print(f"Processing {name_sample}", file=sys.stderr)

                    # Build reference index if not found
                    index_file = Path(dir_ref_pth) / f"{name_fasta}_methylpy_f.4.bt2"
                    print("Index_file: {}".format(index_file))
                    fasta_file = Path(dir_ref_pth) / f"{name_fasta}.fa"
                    print("Ref_fasta_file: {}".format(fasta_file))

                    if os.path.exists(index_file):
                        indx_msg = f"\tIndex file found for {name_fasta}\n\t>>>Skipping creating indices<<<"
                        print(indx_msg, file=fh_log)
                        if verbose:
                            print(indx_msg, file=sys.stderr)
                    else:
                        indx_msg = f"\tIndex file not found for {name_fasta}\n\t>>>Creating indices {name_fasta} <<<"
                        #print(indx_msg, file=fh_log)
                        run_build_reference(name_fasta, fasta_file, dir_ref_pth, num_proc, fh_log, verbose=True)
                        if verbose:
                            print(indx_msg, file=sys.stderr)
                    
                    ## process data for later methylation calling
                    bam_path = Path(dir_out_pth) / f"{name_sample}_processed_reads.bam"
                    print("Bam file: {}".format(bam_path))
                    
                    fastq_path = Path(dir_reads_pth) / f"{name_sample}.fastq.gz"
                    print("fastq file: {}".format(fastq_path))
                    
                    stats_path = Path(dir_out_pth) / f"stats_{name_sample}.txt"
                    print("stat file: {}".format(stats_path))
                    
                    mpy_er_path = Path(dir_out_pth) / f"mpy_er_{name_sample}.txt"
                    print("stat file: {}".format(mpy_er_path))
                    
                    fwd_ref_file = Path(dir_ref_pth) / f"{name_fasta}_methylpy_f"
                    print("fwd_ref_file: {}".format(fwd_ref_file))
                    
                    rv_ref_file = Path(dir_ref_pth) / f"{name_fasta}_methylpy_r"
                    print("rv_ref_file: {}".format(rv_ref_file))
                    
                    ref_fasta_pth = Path(dir_ref_pth) / f"{name_fasta}.fa"
                    print("ref_fasta_pth: {}".format(ref_fasta_pth))

                    if bam_path.exists():
                        log_msg = (
                            f"\tAlignment BAM file found for {name_sample}\n"
                            f"\tSkipping ahead to filtering BAM file..."
                        )
                        print(log_msg, file=fh_log)
                        if verbose:
                            print(log_msg, file=sys.stderr)
                    else:
                        log_msg = (
                            f"\tRunning methylpy for {name_sample} using reference amplicon "
                            f"{name_fasta}.fa"
                        )
                        run_single_end_pipeline(name_sample,bam_path, fastq_path, fwd_ref_file, rv_ref_file,
                            ref_fasta_pth, num_proc, dir_out_pth, mpy_er_path,
                            fh_log, verbose=False)
                        if verbose:
                            print(log_msg, file=sys.stderr)

                        
                    #methylpy proccessing    
                                

                    #func 3
                    # start bam operations and filtering
                    
                    filt_msg=f"\tFiltering highly methylated reads"
                    print(filt_msg, file=fh_log)
                    if verbose:
                        print(filt_msg, file=sys.stderr)
                     
                    filtering_er_path = Path(dir_out_pth) / f"filter_er_{name_sample}.txt"
                    print("stat file: {}".format(filtering_er_path))
                    
                    #only include name wiyhout ".bam"
                    filtered_bam_path = Path(dir_out_pth) / f"{name_sample}_filtered_reads"
                    print("Bam file: {}".format(bam_path))
                    
                    
                    bqf_fn(name_sample,bam_path,filtered_bam_path, fasta_file,siteCH ,filT, filtering_er_path,
                            fh_log, verbose=True)
                    
                    
                    # Remove original BAM and index
                    #remove_bam(fh_log,bam_path)
                    
                    # Call methylation
                    filtered_bam = Path(dir_out_pth) / f"{name_sample}_filtered_reads.bam"
                    
                    cms_er_path = Path(dir_out_pth) / f"cms_er_{name_sample}.txt"
                    print("stat file: {}".format(cms_er_path))
                    
                    cm_log_msg = f"\tCalling methylation of filtered alignment file of {name_sample}"
                    print(cm_log_msg, file=fh_log)
                    if verbose:
                        print(cm_log_msg, file=sys.stderr)
                    
                    #calling the fucntion
                    call_methylation_state(name_sample,filtered_bam,
                            ref_fasta_pth, num_proc, dir_out_pth, cms_er_path,
                            fh_log, verbose=True)
                
                    # Calculate coverage and statistics
                    print("\tCalculating methylation percentage, coverage and statistics", file=fh_log)
                    if verbose:
                        print("\tCalculating methylation percentage, coverage and statistics", file=sys.stderr)
                    
                    cal_cov_stats(fasta_file=ref_fasta_pth,name_sample=name_sample,dir_ref_pth=dir_ref_pth,dir_out_pth=dir_out_pth,fh_log=fh_log)
            

                    
                    
        except FileNotFoundError:
            print(f"[ERROR] Master file not found: {master_file}")
            print(f"[ERROR] Master file not found: {master_file}", file=fh_log)
            sys.exit(1)
        except OSError as e:
            print(f"[ERROR] Could not read {master_file}: {e}")
            print(f"[ERROR] Could not read {master_file}: {e}", file=fh_log)
            sys.exit(1)

        
        #should be change to a seperate option
        #need to allow users to order the samples
        #plotting
        print("\tPlotting output matix results", file=fh_log)
        if verbose:
            print("\tPlotting output matix results", file=sys.stderr)
        
        plot_amplicon(ref_dir=dir_ref_pth,dir_out_pth=dir_out_pth)
        
        # Log completion
        completion_msg = "Analysis pipeline completed successfully!"
        print(completion_msg, file=fh_log)
        if verbose:
            print(completion_msg, file=sys.stderr)
    
        
    # Print final success message (outside all context managers)
    print(f"Results saved to: {dir_out_pth}")        
        
        
        
        
        


if __name__ == "__main__":
    
    start_time = time.time()
    # Missing argument checks
    if not dir_reads:
        print("Error >>>>>>>>>>>> Directory containing sequencing reads not provided <<<<<<<<<<<<<<<<<")
        print_help()
    elif not dir_ref:
        print("Error >>>>>>>>>>>> Directory containing reference files not provided <<<<<<<<<<<<<<<<<")
        print_help()
    elif not master_file:
        print("Error >>>>>>>>>>>> Master file not provided <<<<<<<<<<<<<<<<<")
        print_help()
    
    
    # Path existence checks
    if not os.path.isdir(dir_reads):
        print("Error >>>>>>>>>>>> Reads directory does not exist: {} <<<<<<<<<<<<<<<<<".format(dir_reads))
        sys.exit(1)
    if not os.path.isdir(dir_ref):
        print("Error >>>>>>>>>>>> Reference directory does not exist: {} <<<<<<<<<<<<<<<<<".format(dir_ref))
        sys.exit(1)
    if not os.path.isfile(master_file):
        print("Error >>>>>>>>>>>> Master file does not exist: {} <<<<<<<<<<<<<<<<<".format(master_file))
        sys.exit(1)
    
    
    main()
    print("Analysis pipeline End")

    end_time = time.time()
    elapsed_time = end_time - start_time
    
    print(f"Elapsed time: {elapsed_time:.4f} seconds")









