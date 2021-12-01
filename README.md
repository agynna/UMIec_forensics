


# umierrorcorrect

Pipeline for analyzing barcoded STR markers sequencing data with Unique molecular identifiers (UMI). This pipeline combines correction of sequencing errors using umierrorcorrect with typing of forensic genetic markers (STRs) with FDStools. Output is identical to the output of FDStools.

##Installation
Prior to running one most install [AdaterRemoval](https://github.com/MikkelSchubert/adapterremoval) and [Samtools](https://github.com/samtools/samtools) and make (https://www.gnu.org/software/make/)
The pipeline also makes use of specialzied version of [FLASH](https://academic.oup.com/bioinformatics/article/27/21/2957/217265) made by Jerrythafasta[https://github.com/Jerrythafast]: [FLASH-lowercase-overhang](https://github.com/Jerrythafast/FLASH-lowercase-overhang)
Then umierrorcorrect_forensics can be installed by running:
cd umierrorcorrect_forensics
pip install .
./setup_flash.py

Example command to type STRs from paired ends:

umierrorcorrect_forensics/run_umierrorcorrect_forensics.py -r1 data/15-2800M-10ng_S15_L001_R1_001.fastq.gz \
-r2 data/15-2800M-10ng_S15_L001_R2_001.fastq.gz -p -o results -l data/ULTRALibrary.txt -b data/ULTRALibrary.bed -g data/reference_genome/hg38.fa -i data/ultra.ini

Example command for single ends:

umierrorcorrect_forensics/run_umierrorcorrect_forensics.py -r1 data/15-2800M-10ng_S15_L001_R1_001.fastq.gz \
-o results -l data/ULTRALibrary.txt -b data/ULTRALibrary.bed -g data/reference_genome/hg38.fa -i data/ultra.ini
