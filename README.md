# Per Read Variant Calller (PRVC)

**Description**   
This script, PRVC, performs per read variant calling for  MinION long read sequencing data. PRVC is run from the command line and takes a reference file and one or many associated bam files as required inputs as well as an output file as an optional input. PRVC generates a gzipped csv file with read name, read length, read chromosome, read start position, read end position, number of SNPS, number of deletions, number of insertions, and number of splicing gaps for each read.

# PRVC Grapher

**Description**  
This file creates graphs that visualize the distribution of bases and SNPs for the purpose of aiding in other Pai Lab research.

# File Monitor

**Description**  
This script generates a detailed report on the content of a given directory and identifies the location of specified file types. The report consists of two fiels: one containing the owner, path, and size of all folders and files in a directory (and all its sub directories) and the other is a contains owner, filesize (Gb), file type, file path for every existitng file of specified type in the given directory. This script was used to monitor the  storage of the lab directory on the HPCC. 
