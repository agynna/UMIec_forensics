#!/usr/bin/env python3
import pandas as pd

def stutter_names(allele, minimum_repeat_no=1, kind=-1): 
    '''
    By default returns a list of possible minus one 
    stutter alleles, with at least 'minimum_repeat_no' 
    copies of each repeat left (default 1). Uses 
    bracketed STRNaming sequences,  
    i. e. 'CE23_TAGA[15]CAGA[8]_+1T>C'.
    '''
    if kind not in [-1, 1]: 
        ValueError("Kind must be '-1' or +1. To get minus two stutters, run function twice recursivly. ")

    if len(allele.split("_")) > 2:
        ce_string, seq, flanks = allele.split("_")
    else: 
        ce_string, seq = allele.split("_")

    # Generate new CE number
    new_ce_no = float(ce_string[2:])+kind
    if new_ce_no.is_integer():
        new_ce_string = "CE" + str(int(new_ce_no))
    else: 
        new_ce_string = "CE" + str(round(new_ce_no,1))

    # Unpack bracketed sequence  
    seq_reps = seq.split("]")
    seq_split = []
    for rep in seq_reps: 
        if not rep: 
            continue
        seq_split.append(rep.split("["))
    
    # Create stutter brackets 
    make_stutter_for_rep = []
    stuttered_reps = []
    for rep in seq_split:
        repseq = rep[0]
        repno = int(rep[1])
        if repno > minimum_repeat_no: 
            make_stutter_for_rep.append(True)
            repno = str(repno+kind)
            stuttered_reps.append([repseq, repno])
        else:
            make_stutter_for_rep.append(False)
            stuttered_reps.append(rep)

    # Generate stutter sequences 
    stuttered_seqs = []
    for i in range(len(make_stutter_for_rep)): 
        if make_stutter_for_rep[i]: 
            new_seq = ""
            for j in range(len(make_stutter_for_rep)): 
                if i==j: 
                    new_rep = stuttered_reps[j][0] + "[" + stuttered_reps[j][1] + "]"
                else: 
                    new_rep = seq_split[j][0] + "[" + seq_split[j][1] + "]"
                new_seq = new_seq + new_rep
            stuttered_seqs.append(new_seq)
    
    # Assemble str names 
    stutter_alleles = []
    for seq in stuttered_seqs: 
        if len(allele.split("_")) > 2: 
            new_allele = new_ce_string + "_" + seq + "_" + flanks
        else:
            new_allele = new_ce_string + "_" + seq
        stutter_alleles.append(new_allele)
    return stutter_alleles

def get_stutters(str_names, minimum_repeat_no=1, kind=-1):
    '''
    Returns all possible minus one or plus one stutter allels, 
    with at least 'minimum_repeat_no' copies of each repeat 
    left (default 1). 

    Uses bracketed STRNaming sequences,  i. e. 
    'CE23_TAGA[15]CAGA[8]_+1T>C'. Accepts a single string or a list 
    of strings, and then returns a list of stutter allele strings. 
    Alternativly, a Pandas Series or DataFrame returns a DataFrame 
    with main and corresponding stutter allels in separate columns. 
    '''
    if isinstance(str_names, str): 
        return stutter_names(str_names, minimum_repeat_no, kind)
    elif isinstance(str_names, list):
        stutters = []
        for seq in str_names: 
            stutters.extend(stutter_names(seq, minimum_repeat_no, kind))
        return stutters
    elif isinstance(str_names, pd.Series):
        df_stutters = pd.DataFrame(columns=["main_allele", "stutters"])
        for seq in str_names: 
            s = pd.DataFrame(columns=["main_allele", "stutters"])
            s["stutters"] = stutter_names(seq, minimum_repeat_no, kind)
            s["main_allele"] = seq
            df_stutters = pd.concat([df_stutters, s], axis=0)
        return df_stutters
    elif isinstance(str_names, pd.DataFrame): 
        return get_stutters(str_names.iloc[:,0],minimum_repeat_no, kind)
    else: 
        raise ValueError("Input 'str_names' is of unknown type. (string, list or Pandas Series allowed)")
