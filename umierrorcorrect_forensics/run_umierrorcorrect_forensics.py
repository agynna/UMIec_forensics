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
from umierrorcorrect.filter_bam import filter_bam
from umierrorcorrect.get_consensus_statistics import run_get_consensus_statistics
from run_fdstools import run_fdstools
from convert_fastq2bam import fastq2bam
from convert_bam2fastq import bam2fastq

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

def set_args_preprocessing(args, read_name, output_path, **kwargs):
    tmpdir = kwargs.get('tmpdir', None)
    reads_file = kwargs.get('reads_file', None)
    args.tmpdir = tmpdir
    if reads_file:
        args.read1 = reads_file
    args.output_path = output_path
    read_filename = os.path.basename(args.read1)
    args.sample_name = read_name
    args.mode = 'single'   # We do not use the double end read feature of umierrorcorrect
    args.gziptool = 'gzip'
    args.umi_length = 12
    args.spacer_length = 16
    return args

def set_args_umierrorcorrect(args, read_name, bam_file, bed_file):
    args.bam_file = bam_file
    args.bed_file = bed_file
    args.sample_name = read_name
    args.reference_file = 'data/reference_genome/hg38/hg38.fa'
    args.consensus_frequency_threshold = 0.5
    args.indel_frequency_threshold = 0.6
    args.position_threshold = 20
    args.edit_distance_threshold = 1
    args.regions_from_bed = True
    args.include_singletons = False
    args.remove_large_files = False
    return args

def main():
    # If paired ends, combine reads into single reads using FLASH
    if args.p:
        read_name = common_part_of_read_name(os.path.basename(args.read1),
                                             os.path.basename(args.read2))
        output_path = make_outputdir(args.output_path, read_name)
        tmp_dir = tempfile.mkdtemp()
        merged_reads_file = run_flash(args.read1,
                                        args.read2,
                                        args.num_threads,
                                        tmp_dir,
                                        output_path)
        args_preprocessing = set_args_preprocessing(args, read_name,
                                                    output_path,
                                                    tmpdir=tmp_dir,
                                                    reads_file=merged_reads_file)
    else:
        read_name = os.path.basename(args.read1.split('.',1)[0])
        output_path = make_outputdir(args.output_path, read_name)
        args_preprocessing = set_args_preprocessing(args, read_name,
                                                    output_path)

    # Preprocessing
    (fastq_file, nseqs) = run_preprocessing(args_preprocessing)
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
    bam_file = fastq2bam(output_path, bam_file, bed_file,
                         args.library_file, args.num_threads)

    # Run UMIerrorcorrect
    args_umierrrorcorrect = set_args_umierrorcorrect(args, read_name, bam_file, bed_file)
    run_umi_errorcorrect(args_umierrrorcorrect)

    consensus_reads_file = os.path.join(output_path, read_name + '_consensus_reads.bam')
    filtered_reads_file = os.path.join(output_path, read_name + '_filtered_consensus_reads.bam')
    consensus_cutoff = 20
    filter_bam(consensus_reads_file, filtered_reads_file, consensus_cutoff)
    # Calculate UMIerrorcorrect statistics
    hist_file = os.path.join(output_path, read_name + '.hist')
    run_get_consensus_statistics(output_path,
                                 consensus_reads_file,
                                 hist_file,
                                 samplename=read_name)

    # Convert to Fastq file and run FDStools
    consensus_fastq_file = os.path.join(output_path, read_name + '_filtered_consensus_reads.fq')
    bam2fastq(filtered_reads_file, consensus_fastq_file)
    run_fdstools(consensus_fastq_file, args.library_file, args.ini_file, output_path)


if __name__ == '__main__':
    args = parseArgs()
    main()
