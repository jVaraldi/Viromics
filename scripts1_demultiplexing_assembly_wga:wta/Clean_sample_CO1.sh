SAMPLE=$1
OUTPUT=$2

cd /beegfs/data/varaldi/VIROMICS/data/reads1

# trim reads 1 on the presence of the expected Trans-plex/Omni-plex adaptor at the beginning and remove 3' ends if the RC of  Trans-plex/Omni-plex adaptor is found. Only ^TGTGGTGTGTTGGGTGTGTTTGG is required but untrimmed reads are removed. 2 mismatches per full adaptor is tolerated without gaps. Overlap min = 18bp (out of 23).

#####################################
# trimm primer HCO forward in read 1
#####################################
cutadapt -e 0.09 --no-indels -g "TAAACTTCAGGGTGACCAAAAAATCA;min_overlap=18" --trimmed-only -o  trimmed-$SAMPLE.1bf.fastq  trimmed-$SAMPLE.1.fastq > log_adaptor_R1f_$SAMPLE.out
n1a=$(grep -c ^@MSQ-M01173 trimmed-$SAMPLE.1.fastq)
n1b=$(grep -c ^@MSQ-M01173 trimmed-$SAMPLE.1bf.fastq)

# and trimm on quality . Very stringent to be sure to have full high quality CO1 sequences
java -jar ~/TOOLS/Trimmomatic-0.38/trimmomatic-0.38.jar SE trimmed-$SAMPLE.1bf.fastq trimmed-$SAMPLE.1bf_clean.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:10:15 MINLEN:220
n1c=$(grep -c ^@MSQ-M01173 trimmed-$SAMPLE.1bf_clean.fastq)

#####################################
# trimm primer LCO forward in read 1
#####################################
cutadapt -e 0.09 --no-indels -g "GGTCAACAAATCATAAAGATATTGG;min_overlap=18" --trimmed-only -o  trimmed-$SAMPLE.1br.fastq  trimmed-$SAMPLE.1.fastq > log_adaptor_R1r_$SAMPLE.out
n1d=$(grep -c ^@MSQ-M01173 trimmed-$SAMPLE.1br.fastq)

# and trimm on quality . Very stringent to be sure to have full high quality CO1 sequences
java -jar ~/TOOLS/Trimmomatic-0.38/trimmomatic-0.38.jar SE trimmed-$SAMPLE.1br.fastq trimmed-$SAMPLE.1br_clean.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:10:15 MINLEN:220
n1e=$(grep -c ^@MSQ-M01173 trimmed-$SAMPLE.1br_clean.fastq)


cd /beegfs/data/varaldi/VIROMICS/data/reads2

#####################################
# trimm primer HCO forward in read 2
#####################################

cutadapt -e 0.09 --no-indels -g "TAAACTTCAGGGTGACCAAAAAATCA;min_overlap=18" --trimmed-only -o  trimmed-$SAMPLE.2bf.fastq  trimmed-$SAMPLE.2.fastq > log_adaptor_R2f_$SAMPLE.out
n2a=$(grep -c ^@MSQ-M01173 trimmed-$SAMPLE.2.fastq)
n2b=$(grep -c ^@MSQ-M01173 trimmed-$SAMPLE.2bf.fastq)

# and trimm on quality . Very stringent to be sure to have full high quality CO1 sequences
java -jar ~/TOOLS/Trimmomatic-0.38/trimmomatic-0.38.jar SE trimmed-$SAMPLE.2bf.fastq trimmed-$SAMPLE.2bf_clean.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:10:15 MINLEN:220
n2c=$(grep -c ^@MSQ-M01173 trimmed-$SAMPLE.2bf_clean.fastq)

#####################################
# trimm primer HCO forward in read 2
#####################################

cutadapt -e 0.09 --no-indels -g "GGTCAACAAATCATAAAGATATTGG;min_overlap=18" --trimmed-only -o  trimmed-$SAMPLE.2br.fastq  trimmed-$SAMPLE.2.fastq > log_adaptor_R2r_$SAMPLE.out
n2d=$(grep -c ^@MSQ-M01173 trimmed-$SAMPLE.2br.fastq)

# and trimm on quality . Very stringent to be sure to have full high quality CO1 sequences
java -jar ~/TOOLS/Trimmomatic-0.38/trimmomatic-0.38.jar SE trimmed-$SAMPLE.2br.fastq trimmed-$SAMPLE.2br_clean.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:10:15 MINLEN:220
n2e=$(grep -c ^@MSQ-M01173 trimmed-$SAMPLE.2br_clean.fastq)


##########################################################################
#subset Read 2 to Read 1 set to got paired end and in the same order)
##########################################################################

# https://github.com/enormandeau/Scripts/blob/master/fastqCombinePairedEnd.pys
cd ../

# merge both combination forward = reverse or reverse - forward
python2.7 ~/TOOLS/fastqCombinePairedEnd.py ./reads1/trimmed-$SAMPLE.1bf_clean.fastq ./reads2/trimmed-$SAMPLE.2br_clean.fastq ' '

np_fr=$(grep -c ^@MSQ-M01173 ./reads1/trimmed-$SAMPLE.1bf_clean.fastq_pairs_R1.fastq)
single_fr=$(grep -c ^@MSQ-M01173 ./reads1/trimmed-$SAMPLE.1bf_clean.fastq_singles.fastq)


python2.7 ~/TOOLS/fastqCombinePairedEnd.py ./reads1/trimmed-$SAMPLE.1br_clean.fastq ./reads2/trimmed-$SAMPLE.2bf_clean.fastq ' '

np_rf=$(grep -c ^@MSQ-M01173 ./reads1/trimmed-$SAMPLE.1br_clean.fastq_pairs_R1.fastq)
single_rf=$(grep -c ^@MSQ-M01173 ./reads1/trimmed-$SAMPLE.1br_clean.fastq_singles.fastq)

##########################################################################
#merge Read 1 and Read 2 
##########################################################################
#python2.7 ~/TOOLS/merge_paired_reads.py ./reads1/trimmed-$SAMPLE.1bf_clean.fastq_pairs_R1.fastq ./reads2/trimmed-$SAMPLE.2br_clean.fastq_pairs_R2.fastq

#python2.7 ~/TOOLS/merge_paired_reads.py ./reads1/trimmed-$SAMPLE.1br_clean.fastq_pairs_R1.fastq ./reads2/trimmed-$SAMPLE.2bf_clean.fastq_pairs_R2.fastq

# reverse complement the LCO-HCO merged reads and cat with the HCO-LCO merged reads
#python2.7 ~/TOOLS/reverse.complement.py ./reads1/trimmed-$SAMPLE.1br_clean.fastq_pairs_R1.fastq_merged.fasta ./reads1/trimmed-$SAMPLE.1br_clean.fastq_pairs_R1.fastq_mergedRC.fasta

#cat ./reads1/trimmed-$SAMPLE.1bf_clean.fastq_pairs_R1.fastq_merged.fas ./reads1/trimmed-$SAMPLE.1br_clean.fastq_pairs_R1.fastq_mergedRC.fasta > ./CO1/$SAMPLE_CO1.fa

##########################################################################
# mv to destination directories
##########################################################################


mv ./reads1/trimmed-$SAMPLE.1bf_clean.fastq_pairs_R1.fastq ./pairs/
mv ./reads2/trimmed-$SAMPLE.2bf_clean.fastq_pairs_R2.fastq ./pairs/

mv ./reads1/trimmed-$SAMPLE.1br_clean.fastq_pairs_R1.fastq ./pairs/
mv ./reads2/trimmed-$SAMPLE.2br_clean.fastq_pairs_R2.fastq ./pairs/

mv ./reads1/trimmed-$SAMPLE.1bf_clean.fastq_singles.fastq ./pairs/
mv ./reads1/trimmed-$SAMPLE.1br_clean.fastq_singles.fastq ./pairs/


##########################################################################
# Write summary statistics 
##########################################################################

echo "$SAMPLE $n1a $n1b $n1c $n1d $n1e $n2a $n2b $n2c $n2d $n2e $np_fr $np_rf $single_fr $single_rf" >> $OUTPUT

