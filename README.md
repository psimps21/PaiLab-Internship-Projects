# Per Read Variant Calller (PRVC)

**Description**   
This script, PRVC, performs per read variant calling for  MinION long read sequencing data. PRVC is run from the command line and takes a reference file and one or many associated bam files as required inputs as well as an output file as on optional input. PRVC generates a gzipped csv file with read name, read length, read chromosome, read start position, read end position, number of SNPS, number of deletions, number of insertions, and number of splicing gaps for each read.

# PRVC Grapher

**Description**  
This file creates graphs that visualize the distribution of bases and SNPs for the purpose of aiding other research in Pai Lab.
