cd /beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_reads_on_merged_assembly
while read contig; do
  echo "$contig"
  #samtools depth -a -r $contig sample_B.sorted.bam  > $contig'_depth.txt'
  samtools view -b sample_B.sorted.bam $contig > ./bam/$contig'.bam'
done <contigs_of_interest_wga_wta.txt

#cat contig*depth.txt > contigs_of_interest_wga_wta.depth.txt

# samtools 1.12 was used
