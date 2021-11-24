import pysam
import os
import sys
import logging
import subprocess
import argparse
def parseArgs():
    parser = argparse.ArgumentParser(description="Runs full FDStools pipeline")
    parser.add_argument('-o', '--output_path', dest='output_path',
                        help='Path to the outut file, required', required=True)
    parser.add_argument('-b', '--bam', dest='bam_file',
                        help='Path to the BAM file, required', required=True)
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting convert bam to fastq')
    return(args)

def bam2fastq(infile, outfile):
    with open(outfile, 'w') as of:
        subprocess.run(['samtools', 'fastq',
                        infile],
                        stdout=of,
                        check=True)
    # pysam.fastq(infile,
    #             '-o', outfile,
    #             catch_stdout=False)

def main(args):
    bam2fastq(args.bam_file, args.output_path)

if __name__ == '__main__':
    args = parseArgs()
