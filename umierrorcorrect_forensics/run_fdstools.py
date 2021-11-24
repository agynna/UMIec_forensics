#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import logging
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


def run_fdstools(fastq_file, library_file, ini_file):
    subprocess.run(['fdstools', 'pipeline',
                    ini_file,
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
    run_fdstools(args.fastq_file, args.library_file, args.ini_file)

if __name__ == '__main__':
    args = parseArgs()
    main(args)
