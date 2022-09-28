#!/usr/bin/env python3
import argparse
from pydoc import describe
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
from run_umifilter import run_umifilter
from umierrorcorrect_forensics.tools.uncollapse_reads import uncollapse_reads
from umierrorcorrect_forensics.tools.downsample import downsample_reads


def parseArgs():
    parser = argparse.ArgumentParser(description="TSSV seperation of Simsenseq sequencing of forensics markers ")
    parser.add_argument('-o', '--output_path', dest='output_path',
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-r1', '--read1', dest='read1',
                        help='Path to first FASTQ file, R1, required', required=True)
    parser.add_argument('-r2', '--read2', dest='read2',
                        help='Path to second FASTQ file, R2 if applicable')
    parser.add_argument('-l', '--library', dest='library_file',
                        help='Path to the library file for FDStools, required', required=True)
    parser.add_argument('-b', '--bed', dest='bed_file',
                        help='Path to the marker definition file for TSSV, required', required=True)
    parser.add_argument('-i', '--ini', dest='ini_file',
                        help='Path to first FDStools ini file, required', required=True)
    parser.add_argument('-g', '--reference', dest='reference_file',
                        help='reference genome', required=True)
    parser.add_argument('-p', help='If fastq is paired', action='store_true')                    
    parser.add_argument('-c', '--consensus_method', dest='consensus_method',
                        help="Method for consensus generation. One of 'most_common', 'position' or 'MSA'. \
                            [default = %(default)s]", default="most_common")
    parser.add_argument('-cf', '--consensus_frequency_threshold', dest='consensus_frequency_threshold', 
                        help="Minimum proportion of the majority sequence (or base) for consensus to be called. \
                            Only for 'position' and 'most_common' methods. [default = %(default)s]", 
                            type=float, default=0.5)
    parser.add_argument('-um', '--umi_member_threshold', dest='umi_member_threshold', type=int, 
                        help='Minimum number of members in an UMI family. [default = %(default)s]', default=3)
    parser.add_argument('-fm', '--filter_model', dest='filter_model', 
                        help='Path to model for filtering UMI families. No filtering takes place if unset.')
    parser.add_argument('-fth', '--filter_threshold', dest='filter_threshold', type=float,
                        help='Probability threshold for the filtering model. [default=%(default)s]', default=0.5)
    parser.add_argument('-t', '--num_threads', dest='num_threads',
                        help='Number of threads to run the program on. [default=%(default)s]', default='2')
    parser.add_argument('-u', '--uncollapse', dest='uncollapse', 
                        help='Provide uncollapsed FDStools output, useful for diversity evaluation', action='store_true')
    parser.add_argument('--downsample', dest='downsample', type=float,
                        help='Downsample the number of reads by this fraction. Useful for evaluation.')
    parser.add_argument('--sample_seed', dest='seed', type=int,
                        help='Seed for downsampling. For reproducability.')
    parser.add_argument('--keep', dest='keep_large_files', help='Keep large files.', action='store_true')
    
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    logging.info('Starting UMIerrorcorrect forensics')
    return(args)

def check_paths(paths):
    """Check if input files exists."""
    for p in paths:
        if p:
            if not os.path.isfile(p):
                raise FileNotFoundError(p)

def common_read_name(r1, r2):
    """
    Return the longest common part of read file names. First tries 
    to use the part before the first underscore. If those aren't 
    identical, it will use the longest common part from the start. 
    If the names are all different, return the first file name instead.
    """
    if r1.split("_")[0] == r2.split("_")[0]:
        common_part = r1.split("_")[0]
    else: 
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
    args.sample_name = read_name
    args.mode = 'single'   # We do not use the double end read feature of umierrorcorrect
    args.gziptool = 'gzip'
    args.umi_length = 12
    args.spacer_length = 16
    return args

def set_args_umierrorcorrect(args, read_name, bam_file):
    args.bam_file = bam_file
    args.sample_name = read_name
    args.indel_frequency_threshold = 0.6
    args.position_threshold = 20
    args.edit_distance_threshold = 1
    args.regions_from_bed = True
    args.include_singletons = False
    args.remove_large_files = False
    if args.filter_model:
        args.output_json = True 
        args.consensus_frequency_threshold = 0.1
    else:
        args.output_json = False
    if args.umi_member_threshold < 2:
        args.include_singletons = True
    return args

def main():
    check_paths([args.read1, args.read2, args.bed_file, args.ini_file, args.bed_file, args.filter_model])
    # If asked to downsample reads, do that and save into temp folder. 
    if args.downsample: 
        logging.info('Downsampling input reads by factor ' + str(args.downsample))
        tmp_dir = tempfile.mkdtemp()
        if args.p: 
            read1, read2 = downsample_reads(frac=args.downsample, 
                                            read1=args.read1, 
                                            read2=args.read2, 
                                            output_path=tmp_dir, 
                                            seed = args.seed)
        else: 
            read1 = downsample_reads(frac=args.downsample, 
                                     read1=args.read1, 
                                      output_path=tmp_dir,
                                      seed = args.seed)
        logging.info('Selected reads saved in ' + read1 + ' etc.')
    else: 
        read1 = args.read1
        read2 = args.read2


    # If paired ends, combine reads into single reads using FLASH
    if args.p:
        read_name = common_read_name(os.path.basename(read1),
                                     os.path.basename(read2))
        output_path = make_outputdir(args.output_path, read_name)
        if not(args.downsample):
            tmp_dir = tempfile.mkdtemp()
        merged_reads_file = run_flash(read1,
                                      read2,
                                      args.num_threads,
                                      tmp_dir,
                                      output_path)
        args_preprocessing = set_args_preprocessing(args, read_name,
                                                    output_path,
                                                    tmpdir=tmp_dir,
                                                    reads_file=merged_reads_file)
    else:
        read_name = os.path.basename(read1.split('.',1)[0])
        output_path = make_outputdir(args.output_path, read_name)
        args_preprocessing = set_args_preprocessing(args, read_name,
                                                    output_path)

    # Preprocessing using the UMIec preprocessor 
    (fastq_file, nseqs) = run_preprocessing(args_preprocessing)
    fastq_file = fastq_file[0]
    if args.p or args.downsample:
        if args.keep_large_files:
            shutil.copytree(tmp_dir, os.path.join(output_path, "preprocessing_output"))
        else: 
            shutil.rmtree(tmp_dir)

    # Alignment of reads to markers by TSSV
    plot_qc_stats = False # Requires pandas, seaborn & matplotlib.
    tssv_output_path = os.path.join(output_path, "tssv_output")
    run_tssv(fastq_file, args.library_file, args.num_threads, tssv_output_path, plot_qc_stats)

    # Convert to fastq data to BAM file
    # Paired end reads are already trimmed before FLASH, so skip trimming here in that case.
    bam_file = os.path.join(output_path, read_name + '.bam')
    trim_flanks = not args.p
    bam_file = fastq2bam(tssv_output_path, bam_file, args.bed_file,
                         args.library_file, trim_flanks,
                         args.num_threads)
    if not args.keep_large_files:
        shutil.rmtree(tssv_output_path)

    # Run UMIerrorcorrects
    args_umierrrorcorrect = set_args_umierrorcorrect(args, read_name, bam_file)
    run_umi_errorcorrect(args_umierrrorcorrect)

    json_file_path = os.path.join(output_path, read_name + '_umi_families.json')
    consensus_bam_file = os.path.join(output_path, read_name + '_consensus_reads.bam')
    
    if args.filter_model:
        mlfilter_bam_file = os.path.join(output_path, read_name + '_mlfiltered_reads.bam')
        consensus_bam_file = run_umifilter(consensus_bam_file, json_file_path, 
                                    args.filter_model, mlfilter_bam_file,
                                    args.filter_threshold)

    #TODO: Test to include singletons and set umi_member_cutoff = 1. Build new model. 

    filtered_bam_file = os.path.join(output_path, read_name + '_filtered_consensus_reads.bam')
    filter_bam(consensus_bam_file, filtered_bam_file, args.umi_member_threshold)

    # Calculate UMIerrorcorrect statistics
    stats_file = os.path.join(output_path, read_name + '.hist')
    run_get_consensus_statistics(output_path,
                                 consensus_bam_file,
                                 stats_file,
                                 True,
                                 read_name)

    # Convert to Fastq file and run FDStools
    consensus_fastq_file = os.path.join(output_path, read_name + '_filtered_consensus_reads.fq')
    bam2fastq(filtered_bam_file, consensus_fastq_file, num_threads=args.num_threads)
    run_fdstools(consensus_fastq_file, args.library_file, args.ini_file, output_path, verbose=True)
    logging.info('Finished generating consensus sequences!')

    # If desired, the pipeline can output uncollapsed reads files that can be used for diagnosis. 
    if args.uncollapse: 
        logging.info("Generating uncollapsed read files (-u option activated)...")
        uncollapsed_path = os.path.join(output_path, read_name + '_uncollapsed_consensus_reads.fq')
        uncollapse_reads(consensus_fastq_file, uncollapsed_path)
        run_fdstools(uncollapsed_path, args.library_file, args.ini_file, output_path, verbose=False)
        logging.info("Finished generating uncollapsed read files! ")


if __name__ == '__main__':
    args = parseArgs()
    main()
