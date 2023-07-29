#! /bin/bash 
###############
# for CO1 data
###############
date
echo "start CO1"

#create output summary stats file
log2=/beegfs/data/varaldi/VIROMICS/data/summary_stats_CO1.txt
rm $log2
echo "SAMPLE R1 R1-HCO R1-HCO-clean R1-LCO R1-LCO-clean R2 R2-HCO R2-HCO-clean R2-LCO R2-LCO-clean HCO-LCO LCO-HCO HCO-LCOsingle LCO-HCOsingle" > $log2

while read line 
do 
echo "$line\n" 
sh Clean_sample_CO1.sh $line $log2
done < /beegfs/data/varaldi/VIROMICS/data/sample_names_CO1.txt

date
echo "end"
