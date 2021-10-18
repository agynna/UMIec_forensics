import pysam
import os
import logging
import subprocess
import argparse
def parseArgs():
    parser = argparse.ArgumentParser(description="Runs full FDStools pipeline")
    parser.add_argument('-o', '--output_path', dest='output_path',
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-f', '--fastq', dest='fastq_file',
                        help='Path to first FASTQ file, R1, required', required=True)
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting convert fastq to bam')
    return(args)

if __name__ == '__main__':
    args = parseArgs()
