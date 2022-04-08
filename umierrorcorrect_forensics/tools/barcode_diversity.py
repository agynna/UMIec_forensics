#!/usr/bin/env python3
from unittest.util import _Mismatch
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse
import sys
import os
import re
import math
from fdstools.tools import tssv

def parse_arg():
    parser = argparse.ArgumentParser(description = 'Creates umi diversity plot')
    parser.add_argument("-i", "--infile", dest = "infolder", help = "Folder created by FDStools tssv")
    parser.add_argument("-o", "--outfolder", dest = "outfolder", help = "Path to save diversity plot, optional")
    parser.add_argument("-c", "--csvfile", dest = "csvfile", help = "Path to save csv file, optional")
    parser.add_argument("-l", "--libfile", dest = "libfile", help = "FDStools library file, optional, if given sequences will be trimmed to flanking regions")
    args = parser.parse_args(sys.argv[1:])
    return args

def get_family_sizes(filepath):
    '''
    Reads a paired.fq file and extracts the sizes and sequences 
    of all UMI families. Returns a list of sizes and a DataFrame 
    of sizes and sequences.
    '''
    family_sizes = []
    family_seqs = []
    with open(filepath + "paired.fq") as fh:
        for line in fh:
            if line.startswith("@Consensus"):
                family_sizes.append(int(line.split("=")[-1].rstrip()))
                family_seqs.append(fh.readline()[:-1])
    df_family_sizes = pd.DataFrame({"size": family_sizes, 
                                    "sequence": family_seqs})
                                    
    return family_sizes, df_family_sizes

def revcomp(in_seq):
    pairs = {"A":"T", "C":"G", "G":"C","T":"A"}
    rev_seq = ""
    for index in range(len(in_seq) - 1, 0, -1):
        rev_seq += pairs[in_seq[index]]
    return rev_seq

def make_fdstools_library(library_file): 
    '''
    Reads the flanking regions from a FDStools library file and constructs a 
    minimal TSSV library in its internal representation. 
    '''
    mismatches = 0.1
    lib = dict()
    with open(library_file,"r") as f: 
        l = f.readline()
        while(l):
            if l.startswith("[flanks]"):
                l = f.readline()
                while l:
                    if l.startswith(";"): 
                        pass
                    elif l == "\n": 
                        break
                    elif l.startswith("["):
                        break
                    elif len(re.split("=|,",l)) == 3:
                        l_split = re.split("=|,",l)
                        marker = l_split[0].strip()
                        flank1 = l_split[1].strip()
                        flank2 = revcomp(l_split[2].strip())
                        item = [((flank1, 0), (flank2, 0)), 
                                (math.ceil(len(flank1)*mismatches), math.ceil(len(flank2)*mismatches)), 
                                ""]
                        lib[marker] = item
                    else:
                        print(l)
                        raise Exception("Did not recognize the above line in the library file.")
                    l = f.readline()
            l = f.readline()
    return lib

def trim_sequences(df_longseqs, library_file):
    '''
    Removes the flanks from fastq sequences as specified in library file.  
    This calls the FDStools TSSV internal function to do that. This is a bit 
    hacky, and the FDStools manual explicitly tells you not to do that. 
    Use with caution. 
    Takes DataFrame with column "sequence", returns identical DF with that
    column trimmed. 
    '''
    indel_score = 2
    tssv_library = make_fdstools_library(library_file)
    trimmed_seqs = []
    for seq in df_longseqs["sequence"]:
        list_of_markers_trimmed = tssv.process_sequence(tssv_library, indel_score, False, seq)
        # Extract sequence from tssv return structure
        trimmed_seq = ""
        for record in list_of_markers_trimmed: 
            if record[2]:
                trimmed_seq = record[2]
            elif record[3]:
                trimmed_seq = record[3]
        # For some reason, trimmed seqs come out one nt too long. 
        trimmed_seqs.append(trimmed_seq[0:-1])
    df_out = df_longseqs
    df_out["sequence"] = trimmed_seqs
    return df_out

def create_histo_loop(folder, outfolder=None, csvfile=None):
    folder = folder + "/"
    all_subdirs = [d for d in os.listdir(folder) if os.path.isdir(folder + d)]
    df_families = pd.DataFrame(columns=["marker", "size", "sequence"])

    if outfolder: 
        f, axes = plt.subplots(4, 2, figsize=(15, 15), sharex=False) 
        for marker, ax in zip(all_subdirs, axes.flat):
            family_sizes, df_marker_families = get_family_sizes(f"{folder}/{marker}/")
            print(f"{folder}/{marker}/", ", has ", str(len(family_sizes)), " barcodes. ")
            df_marker_families["marker"] = marker
            if library_file:
                df_marker_families = trim_sequences(df_marker_families, library_file)
            df_families = pd.concat([df_families, df_marker_families])
            
            sns.histplot(family_sizes, ax = ax, bins= 50)
            ax.set_title(marker)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
            ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
            ax.set_xlim(0)
            ax.tick_params(labelrotation=90)
        f.tight_layout()
        f.savefig(outfolder + "/" + "Diversity.png")
    else: 
        for marker in all_subdirs:
            family_sizes, df_marker_families = get_family_sizes(f"{folder}/{marker}/")
            # print(f"{folder}/{marker}/", ", has ", str(len(family_sizes)), " barcodes. ")
            df_marker_families["marker"] = marker
            df_families = pd.concat([df_families, df_marker_families])
    return df_families

def get_barcode_diversity_sample(infolder, outfolder=None, csvfile=None, library_file=None):
    '''
    Return dataframe of family sizes and corresponding sequences for a sample.
    Optionally writes them to a csv file.
    If a library file is give, it will trim the sequences according to 
    the flanks specified there. 
    '''
    df_families = create_histo_loop(infolder, outfolder, csvfile=csvfile)
    if library_file:
        df_families = trim_sequences(df_families, library_file)
    if csvfile: 
        df_families.to_csv(csvfile, index=False)
    return df_families

def get_barcode_diversity_experiment(folder, csvfile=None, library_file=None):
    '''
    Return dataframe of family sizes and corresponding seqences for an experiment. 
    Optionally writes them to a csv file.
    If a library file is give, it will trim the sequences according to 
    the flanks specified there. 
    '''
    df_all_families = pd.DataFrame(columns=["sample", "marker", "size", "sequence"])
    for directory in os.listdir(folder):
        sample_path = os.path.join(folder, directory, "fdstools_results_afterUMIerrorcorrect")
        df_sample_families = get_barcode_diversity_sample(sample_path, csvfile=None, library_file=library_file)
        df_sample_families["sample"] = directory
        df_all_families = (pd.concat([df_all_families, df_sample_families])
                             .astype({"sample": "string",
                                      "marker": "string", 
                                      "size": "int", 
                                      "sequence": "string"}))
    
    if csvfile: 
        df_all_families.to_csv(csvfile, index=False)
    return df_all_families

if __name__ == "__main__":
    '''
    Finds family sizes for UMI super-reads.
    outfolder - if given, saves histograms here.
    csvfile - if given, saves csv file here.
    '''
    args = parse_arg()
    get_barcode_diversity_sample(args.infolder, args.outfolder, args.csvfile, args.libfile)