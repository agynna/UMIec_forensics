import argparse
import subprocess
import sys
import os
from umierrorcorrect_forensics.run_flash import run_flash
from umierrorcorrect.preprocess import run_preprocessing
from umierrorcorrect_forensics.run_tssv import run_tssv
from umierrorcorrect.umi_error_correct import run_umi_errorcorrect
from umierrorcorrect.get_consensus_statistics import run_get_consensus_statistics
from umierrorcorrect_forensics.run_fdstools import run_fdstools

def parseArgs():
    parser = argparse.ArgumentParser(description="TSSV seperation of Simsenseq sequencing of forensics markers ")
    parser.add_argument('-o', '--output_path', dest='output_path',
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-r1', '--read1', dest='read1',
                        help='Path to first FASTQ file, R1, required', required=True)
    parser.add_argument('-r2', '--read2', dest='read2',
                        help='Path to second FASTQ file, R2 if applicable')
    parser.add_argument('-l', '--library', dest='library_file',
                        help='Path to the Library file for TSSV, Required', required=True)
    parser.add_argument('-i', '--ini', dest='ini_file',
                        help='Path to first FDStools ini file, required', required=True)
    parser.add_argument('-t', '--num_threads', dest='num_threads',
                        help='Number of threads to run the program on. Default=%(default)s', default='1')
    parser.add_argument('-p', help='If fastq is paired', action='store_true')
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting UMIerrorcorrect forensics')
    return(args)

if __name__ == '__main__':
    args = parseArgs()
    if args.p:
        run_flash()
    run_preprocessing()
    run_tssv()
    run_umi_errorcorrect()
    run_get_consensus_statistics()
    run_fdstools()
