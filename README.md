# umierrorcorrect

Pipeline for analyzing barcoded STR markers sequencing data with Unique molecular identifiers (UMI). This pipeline combines correction of sequencing errors using umierrorcorrect with typing of forensic genetic markers (STRs) with FDStools. Output is identical to the output of FDStools.

Example command to type STRs from paired ends:
python umierrorcorrect_forensics/run_umierrorcorrect_forensics.py -r1 data/15-2800M-10ng_S15_L001_R1_001.fastq.gz \
-r2 data/15-2800M-10ng_S15_L001_R2_001.fastq.gz -p -o results -l data/ULTRALibrary.txt -i data/ultra.ini

Example command for single ends: 
python umierrorcorrect_forensics/run_umierrorcorrect_forensics.py -r1 data/15-2800M-10ng_S15_L001_R1_001.fastq.gz \
-o results -l data/ULTRALibrary.txt -i data/ultra.ini
