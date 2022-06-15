#!/usr/bin/env python3
import argparse
import random
import sys
import os
import gzip as gz
def parseArgs():
    parser = argparse.ArgumentParser(description="To downsample the number of reads from a SimSenSeq experiment. Useful for evaluation.")
    parser.add_argument('-o', '--output_path', dest='output_path',
                        help='Path to the output directory, otherwise use R1 input path.')
    parser.add_argument('-r1', '--read1', dest='read1',
                        help='Path to first FASTQ file, R1, required', required=True)
    parser.add_argument('-r2', '--read2', dest='read2',
                        help='Path to second FASTQ file, R2 if applicable', default = False)
    parser.add_argument('-f', '--fraction', dest='frac',
                        help='fraction of reads to keep', required=True)
    parser.add_argument('-s', '--seed', dest='seed',
                        help='for reproducability')
    args = parser.parse_args(sys.argv[1:])
    return(args)

def read_fastq(fn):
    """
    A generator that reads gzip compressed FASTA format files and yields each
    block of 4 lines
    
    :param fn: Name of the file with path
    :return: Yields list containing 4 lines, each as str element
    """
    
    block = []
    if fn.endswith(".gz"):
        for line in gz.open(fn):
            block.append(line.decode("utf-8").strip())
            if len(block) == 4:
                yield block
                block = []
    else:
        for line in open(fn):
            block.append(line.strip())
            if len(block) == 4:
                yield block
                block = []

def downsample_reads(frac, read1, read2, output_path, seed):
    frac = float(frac)
    if seed:
        random.seed(int(seed))
    if not(output_path): 
        output_path = os.path.dirname(read1)
    
    if read2:
        file1_name = str.split(os.path.basename(read1),".")[0]
        file2_name = str.split(os.path.basename(read2),".")[0]
        out1_path = os.path.join(output_path,file1_name+"_downsampled.fastq")
        out2_path = os.path.join(output_path,file2_name+"_downsampled.fastq")
        out1 = open( out1_path ,"w")
        out2 = open( out2_path ,"w")
        for r1,r2 in zip(read_fastq(read1),read_fastq(read2)):
            if random.random() < frac:
                for l1,l2 in zip(r1,r2):
                    out1.write(l1 + "\n")
                    out2.write(l2 + "\n")
        return out1_path, out2_path
    else:
        file1_name = str.split(os.path.basename(read1),".")[0]
        out1_path = os.path.join(output_path,file1_name+"_downsampled.fastq")
        out1 = open(out1_path,"w")
        for r1 in read_fastq(read1):
            if random.random() < frac:
                for l1 in r1:
                    out1.write(l1 + "\n")
        return out1_path

def main(): 
    downsample_reads(output_path=args.output_path, 
                      frac = args.frac, 
                      read1 = args.read1, 
                      read2 = args.read2, 
                      seed = args.seed)

if __name__ == '__main__':
    args = parseArgs()
    main()
