cd /beegfs/data/varaldi/MACROGEN/HOST_DB

# L. boulardi
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/PQ/AT/PQAT01/PQAT01.1.fsa_nt.gz

# L. heterotoma
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/RI/CB/RICB01/RICB01.1.fsa_nt.gz


# D. melanogaster
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Drosophila_melanogaster/latest_assembly_versions/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz


# D. sim
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/746/395/GCF_016746395.1_Prin_Dsim_3.0/GCF_016746395.1_Prin_Dsim_3.0_genomic.fna.gz

# D. subobscura
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/008/121/235/GCF_008121235.1_UCBerk_Dsub_1.0/GCF_008121235.1_UCBerk_Dsub_1.0_genomic.fna.gz

# D. hydei
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/285/905/GCF_003285905.1_DhydRS2/GCF_003285905.1_DhydRS2_genomic.fna.gz

# D. immigrans
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/153/375/GCA_018153375.1_ASM1815337v1/GCA_018153375.1_ASM1815337v1_genomic.fna.gz

# D. obscura
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/217/835/GCF_002217835.1_Dobs_1.0/GCF_002217835.1_Dobs_1.0_genomic.fna.gz

# D. suzukii
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/340/165/GCF_013340165.1_LBDM_Dsuz_2.1.pri/GCF_013340165.1_LBDM_Dsuz_2.1.pri_genomic.fna.gz

# merge files and filter on minimim length 5kb (to remove most possible contaminant, see Breitwieser et al.2019) 
# seqfilter from : https://github.com/clwgg/seqfilter
zcat *.gz | gzip > all.gz
~/TOOLS/seqfilter/seqfilter -n -m 5000 -o all_large_scaff_only.fq -i all.gz
gzip all_large_scaff_only.fq > all_large_scaff_only.fq.gz

 


