# Fragilaria radians RNA-seq data analysis

Commands and functions to process the data and generate figures for the manuscript "De novo transcriptome assembly and analysis of the freshwater araphid diatom Fragilaria radians (the former name Synedra acus subsp. radians) isolated from Lake Baikal" by Galachyants et al.

Author: Yuri Galachyants yuragal@gmail.com

For more information, see the paper at: XXX

#### Please note: the code can only be used to see how the data analyses were performed. It should not be treated as a ready-to-use software for transcriptome analysis.

The main reason for the code not being ready-to-use is that it was used in specific cluster environment. You may need to add/edit paths to program executables and change the locations of input and output files/directories according to your needs. Before sequence similarity search the reference database files (such as UniRef, NCBI NR, InterProScan databases) have to be downloaded and prepared. Of course, you will also need to install/compile the required external software packages such as BBMap, DRAP, Trinity, Velvet/Oases assemblers, InterProScan, BLAST, Trinotate, hmmscan etc. Do not forget to ensure the executables are in your $PATH.

### The data processing includes three main parts: 
 1) Prepare and clean the raw data, generate several variants of de novo transcriptome assembies and evaluate them -- performed in bash.
 2) Annotate assembly -- performed in bash.
 3) Analyze DE-transcripts -- mainly performed in R.

Files in the repo are prefixed by numbers denoting their order in the pipeline.

### Raw data at NCBI SRA
    https://www.ncbi.nlm.nih.gov/bioproject/PRJNA484600

### Final data availability
    Fragilaria radians transcriptome assembly: https://www.ncbi.nlm.nih.gov/nuccore/GGVJ00000000.1
    Annotation of Fragilaria radians transcriptome and supplementary data: https://doi.org/10.6084/m9.figshare.7557296.v3
    