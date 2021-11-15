import pysam
import os
import logging
import subprocess
import argparse
import fastq2sam

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

def fastq2bam(infolder, outfile, bed_file, library_file):
    samfile = [outfile + '.sam'] # TODO Change to tempfile
    chromsomes, pos, fq_dirs = get_chr_str(bed_file)
    write_header(samfile, chromsomes)
    loop_fds_result(infolder, samfile, chromsomes, fq_dirs, pos, library_file)

    subprocess.run(['samtools',
                    'view',
                    '--bam', # Output i BAM format
                    samfile, '>', outfile])
                    

def main(args):
    fastq2bam(args.infolder, args.outfile, args.bed_file, args.lib)

if __name__ == '__main__':
    args = parseArgs()
    main(args)
