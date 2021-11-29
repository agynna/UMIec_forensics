#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import logging
import re
from umierrorcorrect.version import __version__

def parseArgs():
    parser = argparse.ArgumentParser(description="Runs full FDStools pipeline")
    parser.add_argument('-f', '--fastq', dest='fastq_file',
                        help='Path to FASTQ file, required', required=True)
    parser.add_argument('-i', '--ini', dest='ini_file',
                        help='Path to first FDStools ini file, required', required=True)
    parser.add_argument('-l', '--library', dest='library_file',
                        help='Path to the Library file for TSSV, Required', required=True)
    parser.add_argument('-o', '--output_path', dest='output_path',
                        help='Path to the output directory')
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting FDStools')
    return(args)

def customize_ini_file(old_ini_file, output_path):
    """
    The ini file contains the path where tssv puts its output. To customize its
    position, we must rewrite the ini file for each sample.
    """
    new_ini_file = os.path.join(output_path, "ultra_custom.ini")
    default_fds_output_path = "fdstools_pipeline_results"
    new_fds_output_path = os.path.join(output_path, 'fdstools_results_afterUMIerrorcorrect')
    with open(old_ini_file, "r") as source:
        lines = source.readlines()
    with open(new_ini_file, "w") as target:
        for line in lines:
            target.write(re.sub(default_fds_output_path,
                                new_fds_output_path,
                                line))
    return new_ini_file

def run_fdstools(fastq_file, library_file, ini_file, output_path):
    new_ini_file = customize_ini_file(ini_file, output_path)
    subprocess.run(['fdstools', 'pipeline',
                    new_ini_file,
                    '-l', library_file,
                    '-s', fastq_file],
                    check=True)
    logging.info('FDStools case-sample pipeline finished. ')

    infile = os.path.splitext(fastq_file)[0] + '.csv'
    outfile = os.path.splitext(fastq_file)[0] + '_stutter.csv'
    subprocess.run(['fdstools', 'stuttermark',
                    '-i', infile,
                    '-o', outfile,
                    '-l', library_file])
    logging.info('Applied stutter thresholds using stuttermark: ' + outfile)
    return None

def main(args):
    run_fdstools(args.fastq_file, args.library_file, args.ini_file, args.output_path)

if __name__ == '__main__':
    args = parseArgs()
    main(args)
