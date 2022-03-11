#!/usr/bin/env python3
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import argparse
import sys
import os

def parse_arg():
    parser = argparse.ArgumentParser(description = 'Creates umi diversity plot')
    parser.add_argument("-i", "--infile", dest = "infolder", help = "Folder created by FDStools tssv")
    parser.add_argument("-o", "--outfolder", dest = "outfolder", help = "Path to save diversity plot, optional")
    parser.add_argument("-c", "--csvfile", dest = "csvfile", help = "Path to save csv file, optional")
    args = parser.parse_args(sys.argv[1:])
    return args

def get_family_sizes(filepath):
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
    if csvfile: 
        df_families.to_csv(csvfile, index=False)
    return df_families

def get_barcode_diversity_sample(infolder, csvfile=None):
    '''
    Return dataframe of family sizes and corresponding sequences for a sample.
    Optionally writes them to a csv file.
    '''
    df_families = create_histo_loop(infolder, outfolder=None, csvfile=csvfile)
    return df_families

def get_barcode_diversity_experiment(folder, csvfile=None):
    '''
    Return dataframe of family sizes and corresponding seqences for an experiment. 
    Optionally writes them to a csv file.
    '''
    df_all_families = pd.DataFrame(columns=["sample", "marker", "size", "sequence"])
    for directory in os.listdir(folder):
        sample_path = os.path.join(folder, directory, "fdstools_results_afterUMIerrorcorrect")
        df_sample_families = get_barcode_diversity_sample(sample_path)
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
    create_histo_loop(args.infolder, args.outfolder, args.csvfile)