#!/usr/bin/env python3
import pysam
import sys
import os
import logging
import argparse
import tempfile
import fastq2sam as f2s

def parseArgs():
    parser = argparse.ArgumentParser(description="Runs full FDStools pipeline")
    parser.add_argument('-o', '--output_file', dest='outfile',
                        help='BAM file to write to, required', required=True)
    parser.add_argument('-f', '--fastq', dest='infolder',
                        help='Path to the TSSV verbose output folder, required', required=True)
    parser.add_argument('-b', '--bed', dest='bed_file',
                        help='Path to BED file with genomic positions of markers, required', required=True)
    parser.add_argument("-l", "--Library", dest = "lib",
                        help = "FDStools library file with marker definitions", required=True)
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting convert fastq to bam')
    return(args)

def fastq2bam(infolder, outfile, bed_file, library_file, trim_flanks=True, num_threads=1):
    [_, samfile] = tempfile.mkstemp(suffix='.sam', text=True)
    chromosomes, pos, fq_dirs = f2s.get_chr_str(bed_file)
    f2s.write_header(samfile, chromosomes)
    f2s.loop_fds_result(infolder, samfile, chromosomes, fq_dirs, pos, library_file, trim_flanks)
    logging.info('Converted fastq to SAM file: ' + infolder + ' to '+ samfile)
    pysam.view('-b', 
               '--output', outfile, 
               samfile, 
               catch_stdout=False)
    logging.info('Compressed SAM to BAM file: ' + outfile)
    os.remove(samfile)

    # outfile_sorted = outfile + '_sorted'
    pysam.sort('-o', outfile,
               '-@', str(num_threads), 
               outfile)
    logging.info('Sorted BAM file: ' + outfile)
    
    pysam.index(outfile)
    logging.info('Indexed BAM file. Fastq to BAM file conversion complete.')
    return outfile

def main(args):
    fastq2bam(args.infolder, args.outfile, args.bed_file, args.lib)

if __name__ == '__main__':
    args = parseArgs()
    main(args)
