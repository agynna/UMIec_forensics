import sys
import argparse

def parse_arg():
    parser = argparse.ArgumentParser(description = 'Uncollapses umibarcodes')
    parser.add_argument("-i", "--infile", dest = "infile", help = "fastq file with umi count header")
    parser.add_argument("-o", "--outfile", dest = "outfile", help = "fastq outfile name")
    args = parser.parse_args(sys.argv[1:])
    return args

def uncollapse_reads(infile, outfile):
    # Writes the number of reads the size of each barcode family
    count = 0
    with open(infile) as fh, open(outfile,"w") as fh_2:
        fastq_lines = []
        for line in fh:
            fastq_lines.append(line)
            if len(fastq_lines) == 4:
                header_uncolapsed = fastq_lines[0].split("=")[0] + "\n"
                numb_reads_uncolapsed = int(fastq_lines[0].split("=")[1])
                fastq_lines[0] = header_uncolapsed
                count += numb_reads_uncolapsed
                fastq_lines_uncolapsed = fastq_lines * numb_reads_uncolapsed
                for line_uncolapsed in fastq_lines_uncolapsed:
                    fh_2.write(line_uncolapsed)
                fastq_lines = []
        print(f"{count} uncollapsed reads")

def main(args):
    uncollapse_reads(args.infile, args.outfile)

if __name__ == '__main__':
    args = parseArgs()
    main(args)
