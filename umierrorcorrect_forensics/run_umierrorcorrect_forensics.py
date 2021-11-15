import argparse
import subprocess
import sys
import os
import tempfile
import shutil
import logging
from run_flash import run_flash
from umierrorcorrect.preprocess import run_preprocessing
from run_tssv import run_tssv
from umierrorcorrect.umi_error_correct import run_umi_errorcorrect
from umierrorcorrect.get_consensus_statistics import run_get_consensus_statistics
from run_fdstools import run_fdstools
from convert_fastq2bam import fastq2bam

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
    parser.add_argument('-g', '--reference-genome', dest='reference_genome',
                        help='Path to human reference genome, required', required=True)
    parser.add_argument('-t', '--num_threads', dest='num_threads',
                        help='Number of threads to run the program on. Default=%(default)s', default='1')
    parser.add_argument('-p', help='If fastq is paired', action='store_true')
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting UMIerrorcorrect forensics')
    return(args)

def common_part_of_read_name(r1, r2):
    """
    Return the longest common part of file names, beginning from the start. If
    they are all different, return the first file name instead.
    """
    i = len(r2)
    while not r1.startswith(r2[0:i]):
        i = i-1
    if i > 0:
        common_part = r1[0:i]
    else:
        common_part = r1
    return common_part

def make_outputdir(output_path, read_name):
    """
    Creates output directory and returns path.
    """
    output_path = os.path.join(output_path, read_name)
    if not os.path.isdir(output_path):
        os.makedirs(output_path)
        logging.info('Created output directory ' + os.path.abspath(output_path))
    else:
        logging.warning('Output directory "' + os.path.abspath(output_path) + '" exists. May overwrite files.')
    return output_path

def main():
    # If paired ends, combine reads into single reads using FLASH
    if args.p:
        read_name = common_part_of_read_name(os.path.basename(args.read1),
                                             os.path.basename(args.read2)
                                            )
        output_path = make_outputdir(args.output_path, read_name)
        tmp_dir = tempfile.mkdtemp()
        merged_reads_file = run_flash(args.read1, args.read2, args.num_threads,
                                        tmp_dir, output_path)
    else:
        read_name = os.path.basename(args.read1)
        output_path = make_outputdir(args.output_path, read_name)

    # Preprocessing
    read_filename = os.path.basename(args.read1)
    args.sample_name = read_filename.split('.',1)[0]
    args.tmpdir = tmp_dir
    args.mode = 'single'   # We do not use the double end read feature of umierrorcorrect
    args.gziptool = 'gzip'
    args.umi_length = 12
    args.spacer_length = 16
    if args.p:
        args.read1 = merged_reads_file
    args.output_path = output_path
    (fastq_file, nseqs) = run_preprocessing(args)
    fastq_file = fastq_file[0]
    if args.p:
        shutil.rmtree(tmp_dir)

    # Alignment of reads to markers by TSSV
    plot_qc_stats = True # Requires pandas, seaborn & matplotlib.
    run_tssv(fastq_file, args.library_file, args.num_threads, output_path, plot_qc_stats)

    # Convert to fastq data to BAM file
    # Assume for now that the BED position file has the same basename and lives
    # in the same dir as the library file.
    bed_file = os.path.splitext(args.library_file)[0] + '.bed'
    bam_file = os.path.join(output_path, read_name + '.bam')
    bam_file = fastq2bam(output_path, bam_file, bed_file, args.library_file)

    run_umi_errorcorrect()
    run_get_consensus_statistics()
    run_fdstools()


if __name__ == '__main__':
    args = parseArgs()
    main()
