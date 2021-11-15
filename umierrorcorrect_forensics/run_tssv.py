#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import logging
import tempfile
import shutil
from umierrorcorrect.version import __version__

def parseArgs():
    parser = argparse.ArgumentParser(description="TSSV separation of Simsenseq sequencing of forensics markers ")
    parser.add_argument('-o', '--output_path', dest='output_path',
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-r1', '--read1', dest='read1',
                        help='Path to first FASTQ file, R1, required', required=True)
#    parser.add_argument('-r2', '--read2', dest='read2',
#                        help='Path to second FASTQ file, R2 if applicable')
    parser.add_argument('-l', '--library', dest='library_file',
                        help='Path to the Library file for TSSV, Required', required=True)
    parser.add_argument('-t', '--num_threads', dest='num_threads',
                        help='Number of threads to run the program on. Default=%(default)s', default='1')
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting TSSV')
    return(args)

def qc_plot_alignment(tssv_folder, out_folder):
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    logging.getLogger('matplotlib.font_manager').disabled=True
    
    stat_df = pd.read_csv(tssv_folder + "/statistics.csv", sep = "\t")
    #flash_stat = pd.read_csv(args.outfolder + "/flash_out.txt", sep = "\t")
    print(stat_df)
    stat_df.set_index("marker",inplace = True)
    marker_balance = stat_df["tPaired"].dropna()

    ax = sns.barplot(x=marker_balance.index, y=marker_balance)
    ax.set(xlabel='Marker', ylabel='Reads')
    fig = ax.get_figure()
    plt.tight_layout()
    fig.savefig(out_folder + "/balance_marker.png")
    #int(flash_stat.iloc[-9][0].split(":")[-1]),
    ax = sns.barplot(x = ["Total reads", "Unrecognised reads", "Reads to markers"],
                     y = [stat_df["unique_seqs"]["total reads"],
                          stat_df["unique_seqs"]["unrecognised reads"],
                          stat_df["tPaired"].sum()]
                    )
    ax.set(xlabel='Marker', ylabel='Reads')
    fig = ax.get_figure()
    plt.tight_layout()
    fig.savefig(out_folder + "/total_reads.png")
    stat_df.to_csv(out_folder + "/statistics_preumi.csv")


def run_tssv(fastq_file, library_file, num_threads, output_path, plot_qc = False):
    # If fastq file is compressed, unzip to temp file.
    if os.path.splitext(fastq_file)[1] == '.gz':
        (_, unzipped_fastq) = tempfile.mkstemp()
        with open(unzipped_fastq, 'w') as uf:
            gzip_run = subprocess.run(['gunzip',
                                       '--to-stdout',
                                       fastq_file],
                                      stdout = uf
                                     )

        if gzip_run.returncode == 0:
            logging.info('Unzipped fastq file')
        else:
            gzip_run.check_returncode()
        fastq_file = unzipped_fastq
        fastq_istempfile = True
    else:
        fastq_istempfile = False

    tssv_run = subprocess.run(['fdstools', 'tssv',
                                '--dir', output_path, # Output dir for verbose output
                                '--indel-score', str(2),
                                '--mismatches', str(0.1),
                                '--num-threads', num_threads,
                                library_file,         # Marker definitions file
                                fastq_file],           # Input file
                                stdout=subprocess.DEVNULL
                              )
    if fastq_istempfile:
        os.remove(fastq_file)
    if tssv_run.returncode == 0:
        logging.info('TSSV finished successfully')
    else:
        tssv_run.check_returncode()

    if plot_qc:
        qc_folder = os.path.join(output_path, 'qc_stats')
        if not os.path.isdir(qc_folder):
            os.mkdir(qc_folder)
        qc_plot_alignment(output_path, qc_folder)
    return None

def main(args):
    run_tssv(args.read1, args.library_file, args.num_threads, args.output_path)
    return None

if __name__ == '__main__':
    args = parseArgs()
    main(args)
