#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import time
import re
from datetime import datetime
from BSAS_runner import run_BSAS

# ===========================================================



def in_arg_parse():
    
    datestring = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    
    parser = argparse.ArgumentParser(
        prog="BSAS_seq",
        description=(
            "Analyze Bisulfite Amplicon sequencing reads - map reads using methylpy, generate allC files,\n"
            " count methylation percentage, coverage,stats and plots data."
        ),
        epilog=(
            "Example usage:\n"
            " python /cluster/pixstor/slotkinr-lab/sandaruwan/seq_pipelines/BSAS/BSA_seq_wrapper/BSA_amp_seq.py\\n"
            "-S /cluster/pixstor/slotkinr-lab/seq_data/BSAS-round20/KS-PHP-amplicon-20\\n"
            "-F /cluster/pixstor/slotkinr-lab/sandaruwan/seq_pipelines/BSAS/BSAS_seq_round20\\n"
            "-M /cluster/pixstor/slotkinr-lab/sandaruwan/seq_pipelines/BSAS/BSAS_seq_round20/PHP_R20_master.txt \\n"
            "-o /cluster/pixstor/slotkinr-lab/sandaruwan/seq_pipelines/BSAS/BSAS_seq_round20/results -v\n\n"    
        ),
        formatter_class=argparse.RawTextHelpFormatter)


    # ====================================== Required arguments ======================================
    parser.add_argument(
        '-S', '--sequences',
        required=True,
        help="Directory containing the sequencing reads files (e.g., /cluster/pixstor/slotkinr-lab/sandaruwan/BSA_round18/KS-SAE-Amplicon-R18).")
    
    parser.add_argument(
        '-F', '--files',
        required=True,
        help="Directory containing the Reference files (e.g., ~/BSA_round18/REF_Files/).")
    
    parser.add_argument(
        '-M', '--master', 
        required=True,
        help="Master file containing the sample names and the reference fasta file names.")
    
    

    # ====================================== optional arguments ======================================

    parser.add_argument(
        '-o', '--output',
        default=f"output_{datestring}",
                        help="Name of the directory to store output files (default: 'output_CURRENTTIME')")
    parser.add_argument('-n', '--numberOfCHSites', type=int, default=5,
                        help="Minimum number of CH sites for mCH level filter (default: 5).")
    parser.add_argument('-t', '--threshold', type=float, default=0.95,
                        help="Threshold of mCH percentage to filter reads (e.g., 96.78 filters >96.78%% CH sites). Use 101 to disable (default: 95).")
    parser.add_argument('-p', '--processors', type=int, default=2,
                        help="Number of processors for methylpy (default: 2).")
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Verbose about parameters and steps (subset of log file).")

    
    return parser.parse_args()


def main_cli():
    
    start_time = time.time()
    args = in_arg_parse()


    # ---- path existence checks ----
    #check input fasta/fastq dir exsits
    if not os.path.isdir(args.sequences):
        sys.exit(f"Error >>> Input files directory does not exist: {args.sequences}")
    #check if reference dir exisits
    if not os.path.isdir(args.files):
        sys.exit(f"Error >>> Reference directory does not exist: {args.ref}")
    #check if master file with amplicn information exisits
    if not os.path.isfile(args.master):
        sys.exit(f"Error >>> bed file of amplicon regions does not exist: {args.region_bed}")
        
    
    print("<<<<<< Starting BSAS amplicon analysis >>>>>>")

    # ---- run pipeline ----
    
    run_BSAS(
        dir_reads = args.sequences,
        dir_ref = args.files,
        dir_out = args.output,
        master_file = args.master,
        filT = args.threshold,
        num_proc = args.processors,
        siteCH = args.numberOfCHSites,
        verbose = args.verbose
    )

    elapsed_time = time.time() - start_time
    
    print("BSAS amplicon analysis pipeline completed")
    print(f"Elapsed time: {elapsed_time:.2f} seconds")


if __name__ == "__main__":
    main_cli()
    
    
    
    

    
    
    
    
    



