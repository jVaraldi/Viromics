# activer conda
#conda activate /beegfs/data/varaldi/myconda

cd /beegfs/data/varaldi/VIROMICS/data

# demultiplex
# trim reads 1 and 2 on the presence of the expected tag id (8bp). 2 differences are allowed.
date
echo "start demultiplexing R1 reads"

cd reads1
cutadapt -e 0.26 --no-indels -g file:../tags_8bp.fasta -o ./trimmed-{name}.1.fastq   ../All_reads1.fastq > log_demultiplex_R1.out

date
echo "start demultiplexing R2 reads"

cd ../reads2
cutadapt -e 0.26 --no-indels -g file:../tags_8bp.fasta -o ./trimmed-{name}.2.fastq   ../All_reads2.fastq > log_demultiplex_R2.out

date
echo "end demultiplexing"
