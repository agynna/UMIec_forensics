import argparse
import subprocess
import sys
import os
import logging
from run_flash import run_flash
from umierrorcorrect.preprocess import run_preprocessing
from run_tssv import run_tssv
from umierrorcorrect.umi_error_correct import run_umi_errorcorrect
from umierrorcorrect.get_consensus_statistics import run_get_consensus_statistics
from run_fdstools import run_fdstools

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

def main(): 
    # Collapse paired end reads into single reads 
    if args.p:
        tmp_path = './tmp/' 
        if not os.path.isdir(tmp_path): 
            os.mkdir(tmp_path) # Move to /temp/ using tempfile? 
            logging.info('Created temporary directory "' + os.path.abspath(tmp_path) + '"')
            
        merged_reads_file = run_flash(args.read1, args.read2, args.num_threads, tmp_path)
        
    # Preprocessing 
    read_filename = os.path.basename(args.read1)
    args.sample_name = read_filename.split('.',1)[0]
    args.tmpdir = tmp_path
    args.mode = 'single'   # We do not use the double end read feature of umierrorcorrect 
    args.gziptool = 'gzip' # Should test if pigz is available, and then serve that instead. 
    args.umi_length = 12
    args.spacer_length = 16
    if args.p:
        args.read1 = merged_reads_file
        
    if not os.path.isdir(args.output_path): 
        os.mkdir(args.output_path) 
        logging.info('Created output directory ' + os.path.abspath(args.output_path)) 
    else: 
        logging.warning('Output directory "' + os.path.abspath(args.output_path) + '" exists. May overwrite files.')
    run_preprocessing(args)
    print(type(args))

    run_tssv()
    
    run_umi_errorcorrect()
    run_get_consensus_statistics()
    run_fdstools()


if __name__ == '__main__':
    args = parseArgs()
    main() 
    
