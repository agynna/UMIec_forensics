#!/usr/bin/env python3
import argparse
import subprocess
import sys
import os
import logging

def parseArgs():
    parser = argparse.ArgumentParser(description="Combines paired reads to one consensus read")
    parser.add_argument('-o', '--output_path', dest='output_path',
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-l', '--log_path', dest='log_path',
                        help='Path to save log file, if different from output')
    parser.add_argument('-r1', '--read1', dest='read1',
                        help='Path to first FASTQ file, R1, required', required=True)
    parser.add_argument('-r2', '--read2', dest='read2',
                        help='Path to second FASTQ file, R2 if applicable')
    parser.add_argument('-t', '--num_threads', dest='num_threads',
                        help='Number of threads to run the program on. Default=%(default)s', default='1')
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting Flash')
    return(args)


def adapter_removal(read1, read2, num_threads, output_path):
    '''
    Removes default Illumina TrueSeq adapters from the reads, which helps
    FLASH with with alignment.
    '''
    trimmed_read1 = os.path.join(output_path, 'trimmed_read1.fastq.gz')
    trimmed_read2 = os.path.join(output_path, 'trimmed_read2.fastq.gz')
    # Adapters from Froste
    adapter1 = 'AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    adapter2 = 'CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCG'
    # Adapters detected
    # adapter1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTNNNNNNNCAATTTCGTATGCCTCTTCTGCTTGA'
    # adapter2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGT'
    subprocess.run(['AdapterRemoval',
                    '--file1', read1,
                    '--file2', read2,
                    '--output1', trimmed_read1,
                    '--output2', trimmed_read2,
                    '--settings', os.path.join(output_path, 'settings.txt'), 
                    '--singleton', os.path.join(output_path, 'single.fasta.gz'), 
                    '--discarded', os.path.join(output_path, 'discarded.fasta.gz'),
                    '--adapter1', adapter1,
                    '--adapter2', adapter2,
                    '--threads', num_threads],
                    check=True)
    return trimmed_read1,trimmed_read2


def run_flash(read1, read2, num_threads, output_path, log_path):
    read1,read2 = adapter_removal(read1, read2, num_threads, output_path)
    read_filename = os.path.basename(read1)
    read_name = read_filename.split('.',1)[0]

    stdout_file = os.path.join(log_path, 'flash_out.txt')
    with open(stdout_file, 'w') as f:
        # Make sure to turn lowercase overhang option -l on!
        flash_run = subprocess.run(['flash',
                                    read1,
                                    read2,
                                    '-t', num_threads,
                                    '-m', str(100),     # Minimum overlap length
                                    '-M', str(300),     # Maximum overlap to be considered in scoring
                                    '-d', output_path,
                                    '-o', read_name,
                                    '-lz'],             # Lowercase overhang + compress output
                                  stdout=f,
                                  stderr=subprocess.STDOUT
                                  )

    if flash_run.returncode == 0:
        logging.info('Flash finished successfully')
    else:
        flash_run.check_returncode()

    output_file = os.path.join(output_path, read_name) + '.extendedFrags.fastq.gz'
    return output_file

def main(args):
    if args.log_path is None:
        log_path = args.output_path
    else:
        log_path = args.log_path
    run_flash(args.read1, args.read2, args.num_threads, log_path)
    return None

if __name__ == '__main__':
    args = parseArgs()
    main(args)
