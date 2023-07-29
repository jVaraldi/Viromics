#!/bin/bash
cd /beegfs/data/varaldi/VIROMICS/data/pairs

#################################
echo "metagenomics data : WGA "
#################################
 
for i in $(seq 57 112);
do
	echo -n "$i \n"
	cat  trimmed-$i'_'*pairs_R1.fastq >> reads1_paired_clean_all_metageno_WGA.fastq
	cat  trimmed-$i'_'*pairs_R2.fastq >> reads2_paired_clean_all_metageno_WGA.fastq

done


# raw data
FILE1=/beegfs/data/varaldi/VIROMICS/data/pairs/reads1_paired_clean_all_metageno_WGA.fastq
FILE2=/beegfs/data/varaldi/VIROMICS/data/pairs/reads2_paired_clean_all_metageno_WGA.fastq
OUTPUT=/beegfs/data/varaldi/VIROMICS/analysis

# check quality of raw reads
fastqc $FILE1 $FILE2 -o $OUTPUT 

#################################
echo "metagenomics data : WTA "
#################################
 
for i in $(seq 113 168);
do
	echo -n "$i \n"
	cat  trimmed-$i'_'*pairs_R1.fastq >> reads1_paired_clean_all_metageno_WTA.fastq
	cat  trimmed-$i'_'*pairs_R2.fastq >> reads2_paired_clean_all_metageno_WTA.fastq

done


# raw data
FILE1=/beegfs/data/varaldi/VIROMICS/data/pairs/reads1_paired_clean_all_metageno_WTA.fastq
FILE2=/beegfs/data/varaldi/VIROMICS/data/pairs/reads2_paired_clean_all_metageno_WTA.fastq
OUTPUT=/beegfs/data/varaldi/VIROMICS/analysis

# check quality of raw reads
fastqc $FILE1 $FILE2 -o $OUTPUT 
