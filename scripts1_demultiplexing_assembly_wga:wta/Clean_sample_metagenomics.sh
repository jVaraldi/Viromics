SAMPLE=$1
SUMMARY_STATS=$2

cd /beegfs/data/varaldi/VIROMICS/data/reads1

#####################################
# trimm adaptor for read 1
#####################################
# trim reads 1 on the presence of the expected Trans-plex/Omni-plex adaptor at the beginning and remove 3' ends if the RC of  Trans-plex/Omni-plex adaptor is found. Only ^TGTGGTGTGTTGGGTGTGTTTGG is required but untrimmed reads are removed. 2 mismatches per full adaptor is tolerated without gaps. Overlap min = 18bp (out of 23).

cutadapt -e 0.09 --no-indels -g "^TGTGGTGTGTTGGGTGTGTTTGG...CCAAACACACCCAACACACCACA;min_overlap=18" --trimmed-only -o  trimmed-$SAMPLE.1b.fastq  trimmed-$SAMPLE.1.fastq > log_adaptor_R1_$SAMPLE.out
# remove remaining reads still containing adaptors
cutadapt -e 0.09 --no-indels -g "TGTGGTGTGTTGGGTGTGTTTGG;min_overlap=8" -a "CCAAACACACCCAACACACCACA;min_overlap=8"  --discard-trimmed -o trimmed-$SAMPLE.1c.fastq trimmed-$SAMPLE.1b.fastq > log_adaptor_R1b_$SAMPLE.out

#########################
# QUALITY TRIMMING for reads 1 
#########################
# from fastqc it appears that 9 bases are GT rich.
#The problem is found on both WGa and WTA datasets. This corresponds to the random primers used in both kits. Maybe they are particlarly GT rich? Anyway, the first 9 bases should probably be removed before final quality trimming and assembly. 
#in addition, we found that the last cycles were of poor quality after 250bp

java -jar ~/TOOLS/Trimmomatic-0.38/trimmomatic-0.38.jar SE trimmed-$SAMPLE.1c.fastq  trimmed-$SAMPLE.1c.clean.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:30 HEADCROP:10 CROP:250


#####################################
# trimm adaptor for read 2
#####################################

cd /beegfs/data/varaldi/VIROMICS/data/reads2

cutadapt -e 0.09 --no-indels -g "^TGTGGTGTGTTGGGTGTGTTTGG...CCAAACACACCCAACACACCACA;min_overlap=18" --trimmed-only -o  trimmed-$SAMPLE.2b.fastq  trimmed-$SAMPLE.2.fastq > log_adaptor_R2_$SAMPLE.out
# remove remaining reads still containing adaptors
cutadapt -e 0.09 --no-indels -g "TGTGGTGTGTTGGGTGTGTTTGG;min_overlap=8" -a "CCAAACACACCCAACACACCACA;min_overlap=8"  --discard-trimmed -o trimmed-$SAMPLE.2c.fastq trimmed-$SAMPLE.2b.fastq > log_adaptor_R2b_$SAMPLE.out

#########################
# QUALITY TRIMMING for reads 2
#########################
# from fastqc it appears that 9 bases are GT rich.
#The problem is found on both WGa and WTA datasets. This corresponds to the random primers used in both kits. Maybe they are particlarly GT rich? Anyway, the first 9 bases should probably be removed before final quality trimming and assembly. 
# in addition, we found that the last cycles were of poor quality after 250bp
java -jar ~/TOOLS/Trimmomatic-0.38/trimmomatic-0.38.jar SE trimmed-$SAMPLE.2c.fastq  trimmed-$SAMPLE.2c.clean.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:30 HEADCROP:10 CROP:250

#########################
# PAIRING
#########################
# https://github.com/enormandeau/Scripts/blob/master/fastqCombinePairedEnd.pys
cd ../
python2.7 ~/TOOLS/fastqCombinePairedEnd.py ./reads1/trimmed-$SAMPLE.1c.clean.fastq ./reads2/trimmed-$SAMPLE.2c.clean.fastq ' '
mv ./reads1/trimmed-$SAMPLE.1c.clean.fastq_pairs_R1.fastq ./pairs/
mv ./reads2/trimmed-$SAMPLE.2c.clean.fastq_pairs_R2.fastq ./pairs/
mv ./reads1/trimmed-$SAMPLE.1c.clean.fastq_singles.fastq ./pairs/


# statistics summary
n1a=$(grep -c ^@MSQ-M01173 reads1/trimmed-$SAMPLE.1.fastq)
n1b=$(grep -c ^@MSQ-M01173 reads2/trimmed-$SAMPLE.2.fastq)

n2a=$(grep -c ^@MSQ-M01173 reads1/trimmed-$SAMPLE.1b.fastq)
n2b=$(grep -c ^@MSQ-M01173 reads2/trimmed-$SAMPLE.2b.fastq)

n3a=$(grep -c ^@MSQ-M01173 reads1/trimmed-$SAMPLE.1c.fastq)
n3b=$(grep -c ^@MSQ-M01173 reads2/trimmed-$SAMPLE.2c.fastq)

n4a=$(grep -c ^@MSQ-M01173 reads1/trimmed-$SAMPLE.1c.clean.fastq)
n4b=$(grep -c ^@MSQ-M01173 reads2/trimmed-$SAMPLE.2c.clean.fastq)

np1=$(grep -c ^@MSQ-M01173 pairs/trimmed-$SAMPLE.1c.clean.fastq_pairs_R1.fastq)
np2=$(grep -c ^@MSQ-M01173 pairs/trimmed-$SAMPLE.2c.clean.fastq_pairs_R2.fastq)
n_single=$(grep -c ^@MSQ-M01173 pairs/trimmed-$SAMPLE.1c.clean.fastq_singles.fastq)



echo "$SAMPLE $n1a $n1b $n2a $n2b $n3a $n3b $n4a $n4b $np1 $np2 $n_single" >> $SUMMARY_STATS

