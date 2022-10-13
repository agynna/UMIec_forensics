import argparse
import sys
import os
import logging
import numpy as np
import pandas as pd
import pickle
import pysam
import lzma

def parseArgs():
    parser = argparse.ArgumentParser(description="Applies an model to filter UMI families.")
    parser.add_argument('-i', '--input_path', dest='input_path',
                        help='Path to the input BAM file with consensus reads.', required=True)
    parser.add_argument('-j', '--json_path', dest='json_path',
                        help='Path to the input json file with UMI family metadata.', required=True)
    parser.add_argument('-o', '--output_path', dest='output_path',
                        help='Path to the output file.', required=True)
    parser.add_argument('-m', '--model_path', dest='model_path',
                        help='Path to a pickle (or .xz) file with a Scikit-learn model to apply.', required=True)
    parser.add_argument('-t', '--threshold', dest='threshold', type=float,
                        help='Classification threshold. [default=%(default)s]', default=0.5)
    args = parser.parse_args(sys.argv[1:])
    logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=logging.DEBUG)
    return(args)

def read_json(json_path):
    """Reads a JSON file with UMI families."""
    df_json = pd.read_json(json_path)
    df_json["UMI"] = df_json["Name"].str.split("_", expand=True)[3]
    df_json["marker"] = df_json["Annotation"].apply(lambda x: x[2])
    # df_json = df_json.rename({"Consensus": "sequence", "Annotation": "marker"}, axis=1).reset_index()
    return df_json

def get_shorter_longer(umi_fam, get_proportions=True): 
    '''Returns number of shorter, same length and longer than consensus members in each UMI family.'''
    n_reads_per_seq = sorted(list(umi_fam["Members"].values()))
    n_pure = n_reads_per_seq[-1]
    if len(n_reads_per_seq) > 1:
        n_unpure = n_reads_per_seq[-2]
    else: 
        n_unpure = 0
    cons_len = len(umi_fam["Consensus"])
    n_variants = len(n_reads_per_seq)

    n_shorter = 0
    n_longer = 0
    n_samelength = 0
    for pair in zip(umi_fam["Members"].keys(), umi_fam["Members"].values()):
        if len(pair[0]) < cons_len: 
            n_shorter = n_shorter + pair[1]
        elif len(pair[0]) > cons_len: 
            n_longer = n_longer + pair[1]
        else: 
            n_samelength = n_samelength + pair[1]
    
    output = [n_pure, n_unpure, cons_len, n_shorter, n_samelength, n_longer, n_variants]
    if get_proportions:
        total_reads = sum(n_reads_per_seq)
        output = [val/total_reads for val in output]
    return output

def get_stutterprops(umi_fam, repeat_length, get_proportions=True): 
    '''Returns proportion of stutters at -1, -2 and +1 positions.'''
    cons_len = len(umi_fam["Consensus"])
    n_minus1 = 0
    n_minus2 = 0
    n_plus1 = 0
    for pair in zip(umi_fam["Members"].keys(), umi_fam["Members"].values()):
        if len(pair[0]) == cons_len-1*repeat_length: 
            n_minus1 = n_minus1 + pair[1]
        elif len(pair[0]) == cons_len-2*repeat_length: 
            n_minus2 = n_minus2 + pair[1]
        elif len(pair[0]) == cons_len+1*repeat_length:
            n_plus1 = n_plus1 + pair[1]

    output = [n_minus1, n_minus2, n_plus1]
    if get_proportions:
        total_reads = sum(list(umi_fam["Members"].values()))
        output = [val/total_reads for val in output]
    return output

def calc_features(df_json):
    """Adds features to the DataFrame with UMI families."""
    df_features = df_json

    # Number of members per UMI group, and number of members compared to average group size
    df_features["total_count"] = df_features["Name"].str.split("=", expand=True)[1].astype(int)
    df_features["normalized_count"] = df_features["total_count"] / df_features["total_count"].mean()

    # Purity, propotion longer, shorter etc
    df_features[["purity", "unpurity", "cons_len", "prop_shorter", "prop_samelength", \
        "prop_longer", "prop_n_variants"]] = df_features.apply(get_shorter_longer, axis=1, result_type="expand")
    df_features["diff_longer_shorter"] = (df_features["prop_longer"]-df_features["prop_shorter"])/df_features["total_count"]

    # Stutter proportion. TODO: Currently only supports markers with repeat length of 4! Fix please. 
    df_features[["prop_minus1", "prop_minus2", "prop_plus1"]] = df_features.apply(get_stutterprops, axis=1, result_type="expand", repeat_length=4)
    return df_features

def apply_filter(df_json, model_path, threshold=0.5):
    """Applies model to UMI families and returns those who passed."""
    df_json = calc_features(df_json)

    # Read ML model. TODO: Change from pickle format to ONNX or PMML. 
    if os.path.splitext(model_path)[1] == ".xz":
        with lzma.open(model_path, "rb") as f: 
            model = pickle.load(f)
    else:
        with open(model_path, "rb") as f: 
            model = pickle.load(f)

    if isinstance(model, dict):
        probabilities = predict_proba_model_dict(model, df_json)
    else:
        probabilities = model.predict_proba(df_json)
    df_out = df_json.loc[probabilities[:,0] >= threshold, :]
    return df_out

def predict_proba_model_dict(model_dict, df):
    """Replaces the predict_proba function when a dict of marker-specific models is used."""
    y = np.empty((len(df), 2), dtype=float)
    marker_list = df.marker.unique()
    for marker in marker_list:
        if not marker in model_dict:
            ValueError(f"Marker name '{marker}' is not present in the user provided ML filter model.")
        y[df["marker"] == marker, :] = model_dict[marker].predict_proba(df.loc[df["marker"] == marker])
    return y

def filter_bamfile(infilename, outfilename, acceptedfams):
    """Reads/writes BAM-file while removing those UMI families that are not in list."""
    n_accepted = 0
    n_dismissed = 0
    with pysam.AlignmentFile(infilename,'rb') as f, pysam.AlignmentFile(outfilename,'wb',template=f) as g:
        reads=f.fetch()
        for read in reads:
            name_contig=read.query_name + "_" + read.reference_name
            if acceptedfams.eq(name_contig).any():
                g.write(read)
                n_accepted += 1
            else:
                n_dismissed += 1
    if n_accepted == 0: 
        logging.warning(f"Zero consensus reads passed ML filter! Instead, {n_dismissed} were dismissed from the consensus BAM file.")
    pysam.index(outfilename)
    logging.info(f"Wrote ML filtered BAM file: {outfilename}. Included {n_accepted}, dismissed {n_dismissed} consensus sequences out of {n_accepted+n_dismissed} in total.")

def run_umifilter(input_path, json_path, model_path, output_path, threshold):
    """Apply ML model to UMI families and write a new model BAM file with passing consensus sequences only."""
    logging.info(f'Applying ML model to filter UMI families. Using model from {model_path} with threshold >= {threshold}.')
    df_json = read_json(json_path)
    df_filtered = apply_filter(df_json, model_path, threshold)
    accepted_UMIfams = df_filtered["Name"].str.cat(df_filtered["Contig"], sep="_")
    filter_bamfile(input_path, output_path, accepted_UMIfams)
    return output_path

def main(args):
    run_umifilter(args.input_path, args.json_path, args.model_path, args.output_path, args.threshold)
    return None

if __name__ == '__main__':
    args = parseArgs()
    main(args)
