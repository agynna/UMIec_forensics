#!/usr/bin/env python3
import argparse
import random
import sys
import gzip as gz
def parseArgs():
    parser = argparse.ArgumentParser(description="TSSV seperation of Simsenseq sequencing of forensics markers ")
    parser.add_argument('-o', '--output_path', dest='output_path',
                        help='Path to the output directory, required', required=True)
    parser.add_argument('-r1', '--read1', dest='read1',
                        help='Path to first FASTQ file, R1, required', required=True)
    parser.add_argument('-r2', '--read2', dest='read2',
                        help='Path to second FASTQ file, R2 if applicable', default = False)
    parser.add_argument('-f', '--fraction', dest='fraction',
                        help='fraction of reads to keep', required=True)
    parser.add_argument('-s', '--seed', dest='seed',
                        help='for reproducibpility', default = 1)
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
def main():
    frac = float(args.fraction)
    random.seed(int(args.seed))
    if args.read2:
         file1_name = args.read1.split(".")[0]
         file2_name = args.read2.split(".")[0]
         out1 = open(file1_name+"_downsampled.fastq","w")
         out2 = open(file2_name+"_downsampled.fastq","w")
         for r1,r2 in zip(read_fastq(args.read1),read_fastq(args.read2)):
             if random.random() < frac:
                 
                 for l1,l2 in zip(r1,r2):
                     out1.write(l1 + "\n")
                     out2.write(l2 + "\n")
    else:
        file1_name = args.read1.split(".")[0]
        out1 = open(file1_name+"_downsampled.fastq","w")
        for r1 in read_fastq(args.read1):
            if random.random() < frac:
                for l1 in r1:
                    out1.write(l1 + "\n")
if __name__ == '__main__':
    args = parseArgs()
    main()
