import argparse
import sys
import re
import jellyfish as jf
def parse_arg():
    parser = argparse.ArgumentParser(description = 'Converts FDStools tssv out put to a single sam format file')
    parser.add_argument("-i", "--infile", dest = "infolder", help = "Folder created by FDStools tssv")
    parser.add_argument("-o", "--outfile", dest = "outfile", help = "SAM file name to write too")
    parser.add_argument("-b", "--bed", dest = "bedfile", help = "BED file to get chromesome and strname from")
    parser.add_argument("-l", "--library", dest = "lib", help = "BED file to get chromesome and strname from")
    args = parser.parse_args(sys.argv[1:])
    return args
def rev_comp(seq):
    """
    Returns reverser compliment of DNA string
    """
    complement = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    seq = seq[::-1]
    return "".join(complement[letter] for letter in seq)
def get_chr_str(bedfile):
    """
    Get names and postions of markers in BED file
    """
    with open(bedfile) as fh:
        fq_dirs = []
        chromsomes = []
        pos = []
        for line in fh:
            line = line.rstrip().split("\t")
            chromsomes.append(line[0])
            pos.append(line[1])
            fq_dirs.append(line[3])
    return (chromsomes,pos,fq_dirs)

def levenshtein(s1, s2):
    """
    Levenshtein distance, implementation from the Algorithm implementation
    wikibook: en.m.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
    This avoids the jellyfish dependency, but adds ~50 s of computation time.
    """
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

def search_for_flank_trim(seq, flank, rev_flank,flank_len,flank_rev_len):
    found = False
    """
    Searches read for either the reverse primer or the forward primer reversed
    """
    for m in re.finditer(flank, seq):

        #print('flank found', m.start(), m.end())
        found = True
        return m.end()
    for m in re.finditer(rev_flank, seq):

        #print('rev flank found', m.start(), m.end())
        found = True
        return m.end()
    if not found:
        for base in range(-1,-len(seq)-flank_len,-1):

            subseq = seq[base-flank_len : base]
            subseq_rev = seq[base-flank_rev_len : base]
            #ham = levenshtein(subseq,flank)
            #ham_rev = levenshtein(subseq_rev,rev_flank)
            ham = jf.levenshtein_distance(subseq,flank)
            ham_rev = jf.levenshtein_distance(subseq_rev,rev_flank)
            if ham <= 2:
                #print('flank found ham',base,base - flank_len, base + flank_len , ham)
                return base
                break
            elif ham_rev <= 2:
                #print('rev flank found ham',base,base - flank_len, base + flank_len , ham)
                return base
                break

def get_flanks_lib(filename):

    """
    Gets flanking sequences from FDStools library file
    """
    flanks = {}
    with open(filename) as fh:
        found_flanks = False
        for line in fh:
            line = line.rstrip()
            if line.startswith("[flanks]"):
                found_flanks = True
                continue
            elif line.startswith("[prefix]"):
                return flanks
            elif line.startswith(";") or line.startswith("_") or line == "":
                continue
            elif found_flanks:
                line = line.replace(" ","").split("=")
                #print(line[0],line[1])
                flanks[line[0]] = line[1].split(",")
def write_header(outfile,chromsomes):
### Write header of sam file
    with open(outfile,"w") as outfh:
        dup = []
        for chrom in chromsomes:
            if chrom not in dup:
                sam_header = f"@SQ	SN:{chrom}	LN:10\n"
                outfh.write(sam_header)
                dup.append(chrom)

def loop_fds_result(indir, outfile, chromsomes, fq_dir, pos, lib):
### Loop through 7 hardcoded different fastq files output of tssv and writes them into sam file
    flanks = get_flanks_lib(lib)
    for dire, chromsome, pos in zip(fq_dir, chromsomes, pos):
        fq2sam(f"{indir}/{dire}/paired.fq", outfile, dire, chromsome, pos, flanks[dire][1], rev_comp(flanks[dire][0]))

def fq2sam(infile, outfile, dir, chrom, pos, flank, flank_rev):
### Write fastq seq and phred score to sam file with hardcoded chrom pos and perfect map qual and trims read att positon from search flank
    with open(infile, "r") as infh, open(outfile, "a") as outfh:
        fq_lines = []
        count = 0
        count_2 = 0
        flank_len = len(flank)
        flank_rev_len = len(flank_rev)
        for line in infh:
            fq_lines.append(line.rstrip().split()[0])
            if len(fq_lines) == 4:
                count_2 += 1
                end = search_for_flank_trim(fq_lines[1], flank, flank_rev, flank_len, flank_rev_len)
                #print("Length of read before triming ", len(fq_lines[1]))
                #print("Flank found at:", end)
                if end is None:
                    fq_lines = []
                    continue
                else:
                    fq_lines[1] = fq_lines[1][:end]
                    fq_lines[3] = fq_lines[3][:end]
                    cigar = len(fq_lines[1])
                    pos_2 = int(pos) + int(cigar)
                    #print("Length of read after triming ", len(fq_lines[1]))
                    sam_line = f"{fq_lines[0][1:]}\t0\t{chrom}\t{pos}\t255\t{cigar}M\t*\t0\t0\t{fq_lines[1]}\t{fq_lines[3]}\tUG:Z:{dir}\n"
                    outfh.write(sam_line)
                    fq_lines = []
                    count += 1

        #print(flank_seq)
        print(f"found {count} reads of {count_2} for {dir}")

def main(args):
   chromsomes, pos, fq_dirs = get_chr_str(args.bedfile)
   write_header(args.outfile, chromsomes)
   loop_fds_result(args.infolder, args.outfile,chromsomes, fq_dirs, pos, args.lib)

if __name__ == "__main__":
    args = parse_arg()
    main(args)
