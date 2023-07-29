# activer conda
conda activate /beegfs/data/varaldi/myconda

# raw data
FILE1a=/beegfs/data/varaldi/VIROMICS/data/Varaldi_Viromics_CGATGT_L001_R1_001.fastq
FILE1b=/beegfs/data/varaldi/VIROMICS/data/Varaldi_Viromics_CGATGT_L001_R2_001.fastq
FILE2a=/beegfs/data/varaldi/VIROMICS/data/Varaldi_Viromics_CGATGT_L001_R1_002.fastq
FILE2b=/beegfs/data/varaldi/VIROMICS/data/Varaldi_Viromics_CGATGT_L001_R2_002.fastq

OUTPUT=/beegfs/data/varaldi/VIROMICS/analysis

# check quality of raw reads

fastqc $FILE1a $FILE1b $FILE2a $FILE2b -o $OUTPUT > log.fastqc




# deactiver conda
conda deactivate

