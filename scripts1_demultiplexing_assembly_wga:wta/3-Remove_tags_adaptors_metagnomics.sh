# apply the cleaning script to remove the identifier and the adaptors
# for each demultiplexed sample

 #! /bin/bash 
###############
# for metagenomics samples
###############
date
echo "start metagenomics"
log1=/beegfs/data/varaldi/VIROMICS/data/summary_stats_metagenomics.txt
#create output summary stats file
rm $log1
echo "samp reads1 reads2 reads1b reads2b reads1c reads2c reads1d reads2c pairedR1 pairedR2 single" >> $log1

while read line 
do 
echo "$line\n" 
sh Clean_sample_metagenomics.sh $line $log1
done < /beegfs/data/varaldi/VIROMICS/data/sample_names_metageno.txt
