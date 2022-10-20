
![UMIerrorcorrect Forensics Logo](https://github.com/FrosteS/umierrorcorrect_forensics/blob/main/Umierrorcorrect_forensics.PNG)

# UMIerrorcorrect Forensics

UMIerrorcorrect Forensics is a pipeline for analyzing barcoded STR markers sequencing data with unique molecular identifiers (UMIs). This pipeline combines correction of sequencing errors using Umierrorcorrect with typing of forensic genetic markers (STRs) with FDStools. Output is identical to the output of FDStools.

## Installation
Prior to running one must install [AdapterRemoval](https://github.com/MikkelSchubert/adapterremoval) and [make](https://www.gnu.org/software/make/). The '-c msa' option (consensus using multiple sequence alignment) requires the [MAFFT](https://mafft.cbrc.jp/alignment/software/) alignment engine. The pipeline also includes a version of [FLASH](https://academic.oup.com/bioinformatics/article/27/21/2957/217265) made by [Jerrythafasta](https://github.com/Jerrythafast): [FLASH-lowercase-overhang](https://github.com/Jerrythafast/FLASH-lowercase-overhang). 

Then umierrorcorrect_forensics can be installed by running:\
```
cd umierrorcorrect_forensics\
pip install .\
./setup_flash.py\
```
## Usages
Example command for single ends:
```
run_umierrorcorrect_forensics.py -r1 data/example-2800M-1ng_R1.fastq.gz \
-o results -l data/ultra_library.txt -b data/ultra_markers.bed -g data/mini_hg38.fa -i data/ultra.ini
```
Example command to type STRs from paired ends:
```
run_umierrorcorrect_forensics.py -r1 data/example-2800M-1ng_R1.fastq.gz \
-r2 data/example-2800M-10ng_R2.fastq.gz -p -o results -l data/ultra_library.txt -b data/ultra_markers.bed -g data/mini_hg38.fa -i data/ultra.ini
```
Example command when using ML filter:
```
run_umierrorcorrect_forensics.py -r1 data/example-2800M-1ng_R1.fastq.gz \
-o results -l data/ultra_library.txt -b data/ultra_markers.bed -g data/mini_hg38.fa -i data/ultra.ini \ 
--consensus_frequency_threshold 0 --umi_member_threshold 2 \
--filter_model data/221012-RFModelDict-separate.xz --filter_threshold 0.95 
```
## Settings
The --library, --bed, --reference and --ini files are mandatory. 

The library, bed and reference genome files should be adapted to the loci included in your assay. The included files describe the original 7-loci Ultra assay. The **library file** is an FDStools library file. It contains settings for each of the loci, and should be self-explanatory. The **bed file** specifies the chromosome and genomic coordinates for each locus. To decrease the size of this package, the included **reference genome** is a minimal version of the human genome including only the regions around the loci of interest. You can either use a full reference genome or create a similar minimal file. The coordinates in the bed file should be adjusted accordingly. 

The **ini file** is an FDStools configuration file. It can likely be used unchanged between different assays, but you might want to change the visualisation settings ("vis") to adapt the html output file to your liking. 

### Consensus method
The consensus method determines how the consensus sequence for each UMI family is generated. The consensus generation is implemented in [Umierrorcorrect](https://github.com/stahlberggroup/umierrorcorrect). 
* *most_common* takes the single most common sequence in each group as consensus. 
* *position* aligns each member sequence from the start, and takes the most common base at each position as the consensus at that position. This method is less well suited when the sequence length varies within the UMI groups, and can then occasionally give consensus sequences that do not occour at all between the members (*purity* = 0). 
* *MSA* uses multiple sequence alignment with MAFFT, using default settings. This is well adapted for families with varying lengths, but is computationally expensive. 

### Filtering
The generated consensus sequences are by default filtered so their UMI families have at least 3 members, and at least 50% of the members support the consensus sequencs. These can be adjusted with the --umi_member_threshold and --consensus_frequency_threshold options. For *most common* and *MSA*, 50% of the members must be fully identical, while for *position* this is calculated at each base position and must be fulfilled at all of them. 

UMIerrorcorrect Forensics also supports applying a machine learning model to filter UMI families dependent on information about the family members. When a ML model is applied, the pipeline calulates various statistics on each family which are used as input to the model. For each UMI family, the model returns a score (probability) between 0 and 1, indicating how reliable its consensus sequence is. Consensus sequences with a score under the given threshold are discarded. 

The model is given as a path to a pickle file, optionally lzma compressed (.xz file extension). The pickle file can be organised in three different ways: 
* A single Scikit-learn model, which is then applied to all UMI families. 
* A dict with the locus name as key, and as values one separate Scikit-learn model for each locus. 
* A dict with the locus name as key, and a list as values. The first position of each list should be a Scikit-learn model, which is first applied to the UMI families of each locus. The second and third positions are then used as parameters in a linear model to adjust the model scores depending on the length of the consensus sequence, according to *-(k\*length+m)*, where the second value is *k* and the third is *m*. This can be used to improve intralocus balance in case the model is biased to shorter (or longer) alleles.
* It is allowed to combine option 2 and 3 at different loci. 

The model itself must be a Scikit-learn pipeline or compatible object, which accepts a Pandas dataframe (through Sklearn-pandas), and outputs a two-class classification score. The probability for a correct sequence should be in the first output column (*i. e.* it should be trained with erroreous sequences having the positive label). This software has been tested with an imbalanced-learn pipeline wrapping a RandomForestClassifier, and a with a CalibratedClassifierCV, but may also work with other types of Scikit-learn models. 

These statistics (*i. e. model features*) are calculated for each UMI family and presented to the model: 
* *marker* - which marker/locus the members align to. 
* *total_count* - the number of members in the family, i. e. the number of reads with the same UMI code
* *normalized_count* - *total_count* divided by the average *total_count* for all UMIs in the same sample 
* *purity* - proportion of members identical to the most common sequence in the family. 
* *unpurity* - proportion of members identical to the second most common sequence in the family. 
* *prop_shorter* - proportion of members shorter than the consensus 
* *prop_samelength* - proportion of members with the same length as the consensus 
* *prop_longer* - proportion of members longer than the consensus
* *prop_minus1* - proportion of members that are one repeat unit shorter than the consensus
* *prop_minus2* - proportion of members that are two repeat units shorter than the consensus
* *prop_plus1* - proportion of members that are one repeat unit longer than the consensus
    * **NOTE:** Currently, only repeat lengths of 4 bases are supported. 
* *diff_longer_shorter* - *prop_longer* minus *prop_shorter*
* *prop_n_variants* - the number of different sequence variants in the family, divided by *total_count*
* *cons_len* the length of the consensus sequence in base pairs
    * **NOTE:** The ML models included in the package do no use the *cons_len* feature. It is only used for the linear adjustment of ML scores that intends to compensate for any length-dependent bias in the ML model itself. 

## Acknowledgements 
Umierrorcorrect Forensics depends on the [Umierrorcorrect](https://github.com/stahlberggroup/umierrorcorrect) software for analysis of UMI tagged sequencing data, developed by Tobias Österlund at the University of Gothenburg, and the [FDStools](https://fdstools.nl/) suite for analysis of forensic STR sequencing data, developed by Jerry Hoogenboom at the Netherlands Forensic Institute. 

UMIerrorcorrect Forensics was developed at the National Forensic Center (NFC) of the Swedish Police Agency in 2021-2022. Please see the included LICENSE.txt file. 