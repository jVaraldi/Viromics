import os 

samples = ["B","C"]
#configfile: "path/to/config.yaml" next time prepare a config file to define sample names list and their path
workdir: "/beegfs/data/varaldi/MACROGEN/scripts"
n_B = list(range(1,521)) # n files after spliting the assembly file (sample B - metagenomic/transcripto mix without wga/wta amplification (all samples except #43)). input for the blastx (mmseqs)
n_C = list(range(1,397)) # n files after spliting the assembly file (sample C - sample #43 only for LhFV sequencing). input for the blastx (mmseqs)
#n_split= list(range(0,78)) # nfiles after spliting the merged assembly
n_split= ['0', '77'] # nfiles after spliting the merged assembly for dag only
n_hisat_build= list(range(1,9)) # n files created by hisat2 build
n_mapping_files= ['57', '168'] # n wga/wta files mapped against the reference panel for dag only 
#n_mapping_files= list(range(57,169)) # n wga/wta files mapped against the reference panel 

rule all:    
    input:
       # "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.1.ht2",
       # "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.2.ht2",
        #"/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.3.ht2",
        #"/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.4.ht2",
        #"/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.5.ht2",
        #"/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.cfq.6.ht2",
        #"/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.7.ht2",
        #"/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.8.ht2",
        #expand("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/{sample}", sample=samples),
        #"/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling",
        #"/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/B_1_dedup_host_removed_on_assembly.sorted_COV_TABLE",
        #"/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_1_dedup_host_removed_on_assembly.sorted_COV_TABLE",
       # "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B_subsampling",
        #expand("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/blastx/file_{sample}.blastx", sample=n_B),
        #expand("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/blastx/file_{sample}.blastx", sample=n_C)
        #"/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling/assembly.blastx"
        #expand("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling/assembly.{sample}.ht2", sample=n_hisat_build)
        expand("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/blastx/assembly.fasta-split-{sample}.blastx", sample=n_split),
        expand("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_wga_wta_on_merged_assembly/sample_{number}.table.txt", number=n_mapping_files),
        expand("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_reads_on_merged_assembly/sample_{sample}.table.txt", sample=samples)
        

rule remove_duplicates:
    input:
        r1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_1.fastq.gz",
        r2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_2.fastq.gz" 
    output:
        r1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_1_dedup.fastq.gz",
        r2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_2_dedup.fastq.gz"
    message: "removing duplicates from sample {input}."
    params:
        pipe_R1="{sample}_1.fastq.gz.pipe",
        pipe_R2="{sample}_2.fastq.gz.pipe",
        pipe_file="{sample}_pipelist.txt"                
    resources: time_min=360, mem_mb=36000, cpus=4 # peaks at 102Gb (see benchmark file)! 
                                                  #next time increase RAM
    benchmark:
        "benchmarks/remove_duplicates/{sample}.tsv"
    conda:
        "envs/fastuniq.yaml"
    shell:
        """
        mkfifo {params.pipe_R1}
        mkfifo {params.pipe_R2}
        gunzip -c {input.r1} > {params.pipe_R1} &
        gunzip -c {input.r2} > {params.pipe_R2} &
        echo '{params.pipe_R1}' > {params.pipe_file}
        echo '{params.pipe_R2}' >> {params.pipe_file}
        fastuniq -i {params.pipe_file} -o {input.r1}.dedup -p {input.r2}.dedup
        gzip {input.r1}.dedup
        gzip {input.r2}.dedup
	mv {input.r1}.dedup.gz {output.r1}
	mv {input.r2}.dedup.gz {output.r2}
        # clean
        rm {params.pipe_R1} {params.pipe_R2} {params.pipe_file}
        """


rule hisat2_index:
    input:
        "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.gz"
    output:
        output1="/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.1.ht2",
        output2="/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.2.ht2",
        output3="/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.3.ht2",
        output4="/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.4.ht2",
        output5="/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.5.ht2",
        output6="/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.6.ht2",
        output7="/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.7.ht2",
        output8="/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.8.ht2"
    message: "constructing host sequence database for mapping (to eliminate them)."
    resources: time_min=180, mem_mb=8000, cpus=4
    benchmark:
        "benchmarks/hisat2_index/tab.tsv"
    conda:
        "envs/hisat2.yaml"
    shell:
        """
        gunzip {input}
        hisat2-build -p 4 /beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq /beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq
        gzip /beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq
        """



rule hisat2_align_remove_host:
    input:
        "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.1.ht2",
        "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.2.ht2",
        "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.3.ht2",
        "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.4.ht2",
        "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.5.ht2",
        "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.6.ht2",
        "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.7.ht2",
        "/beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq.8.ht2",
        reads1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_1_dedup.fastq.gz",
        reads2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_2_dedup.fastq.gz"
    message: "mapping reads on host database from files {input}. Only unmapped reads are outputed."
    resources: time_min=480, mem_mb=64000, cpus=16
    benchmark:
        "benchmarks/hisat2_align_remove_host/{sample}.tsv"
    output:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_1_dedup.fastq.gz_1_host_removed.fastq.gz",
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_1_dedup.fastq.gz_2_host_removed.fastq.gz"
    params: 
        db="DB_contigs_wga_wta/DNA_viruses_wga"
    conda:
        "envs/hisat2.yaml"
    log: "logs/hisat2_mapping_{sample}.log"
    shell:# write only paired-end reads that fail to align concordantly
        'hisat2 --summary-file {input.reads1}.hisat.log -x /beegfs/data/varaldi/MACROGEN/HOST_DB/all_large_scaff_only.fq -1 {input.reads1} -2 {input.reads2} -k 1 -p 16  --un-conc-gz {input.reads1}_%_host_removed.fastq.gz 2> {log}'


rule trimmomatic_pe:
    input:
        r1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_1_dedup.fastq.gz_1_host_removed.fastq.gz",
        r2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_1_dedup.fastq.gz_1_host_removed.fastq.gz"
    output:
        r1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/{sample}_1_dedup_host_removed.fastq.gz",
        r2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/{sample}_2_dedup_host_removed.fastq.gz"
        # reads where trimming entirely removed the mate
        #r1_unpaired="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/{sample}_1_dedup_host_removed.unpaired.fastq.gz",
        #r2_unpaired="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/{sample}_2_dedup_host_removed.unpaired.fastq.gz"
    message: "trimming os sample {input}"
    benchmark:
        "benchmarks/trimmomatic_pe/{sample}.tsv"
    resources: time_min=240, mem_mb=16000, cpus=8
    conda:
        "envs/trimmomatics.yaml"
    log:
        "logs/trimmomatic_{sample}.log"
    threads:
        8
    shell:
        """
        java -jar ~/TOOLS/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads {threads} -trimlog {log} {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} LEADING:10 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:36

        """


rule megahit_assemble:
    input:
        r1 = "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/{sample}_1_dedup_host_removed.fastq.gz",
        r2 = "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/{sample}_2_dedup_host_removed.fastq.gz"
    output:
        #directory(os.path.join("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/", '{sample}','/final.contigs.fa'))
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/{sample}/final.contigs.fa"
    message: "assembling sample {input} (whole dataset)"
    benchmark:
        "benchmarks/megahit_assemble/{sample}.tsv"
    resources: time_min=1200, mem_mb=32000, cpus=8
    conda:
        "envs/megahit.yaml"
    log: "logs/megahit_full_{sample}.log"
    shell:
        'megahit -t 8 -1 {input.r1} -2 {input.r2} -o {output} 2> {log}'
	
rule subsample_reads_C:
    input:
        r1 = "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_1_dedup_host_removed.fastq.gz",
        r2 = "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_2_dedup_host_removed.fastq.gz"
    output:
        r1 = "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_1_dedup_host_removed_subset.fastq.gz",
        r2 = "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_2_dedup_host_removed_subset.fastq.gz"
    message: "sampling reads from sample {input}"
    resources: time_min=120, mem_mb=8000, cpus=4
    benchmark:
        "benchmarks/subsample_reads_C/tab.tsv"
    conda:
        "envs/seqtk.yaml"
    log: "logs/subsample_C.log"
    shell:
        """
        zcat {input.r1} | seqtk sample - -s100 8555009 | gzip - > {output.r1} 2> {log}
        zcat {input.r2} | seqtk sample - -s100 8555009 | gzip - > {output.r2} 2> {log}
        """


rule megahit_assemble_subsample_C:
    input:
        r1 = "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_1_dedup_host_removed_subset.fastq.gz",
        r2 = "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_2_dedup_host_removed_subset.fastq.gz"
    output:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa"
    message: "assembling sample {input} (reduced dataset)"
    benchmark:
        "benchmarks/megahit_assemble_subsample_C/tab.tsv"
    resources: time_min=1200, mem_mb=32000, cpus=8
    conda:
        "envs/megahit.yaml"
    log: "logs/megahit_subsample_C.log"
    shell:
        'megahit -t 8 -1 {input.r1} -2 {input.r2} -o {output} 2> {log}'



rule hisat2_index_assembly_B:
    input:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa"
    output:
        output1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.1.ht2",
        output2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.2.ht2",
        output3="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.3.ht2",
        output4="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.4.ht2",
        output5="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.5.ht2",
        output6="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.6.ht2",
        output7="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.7.ht2",
        output8="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.8.ht2"
    message: "constructing assembly database for mapping (to calculate coverage). sample B."
    resources: time_min=180, mem_mb=8000, cpus=4
    benchmark:
        "benchmarks/hisat2_index/tab.tsv"
    conda:
        "envs/hisat2.yaml"
    shell:
        """
        mkdir -p /beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B
        hisat2-build -p 4 {input} {input}
        """


rule hisat2_index_assembly_C:
    input:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa"
    output:
        output1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.1.ht2",
        output2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.2.ht2",
        output3="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.3.ht2",
        output4="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.4.ht2",
        output5="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.5.ht2",
        output6="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.6.ht2",
        output7="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.7.ht2",
        output8="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.8.ht2"
    message: "constructing assembly database for mapping (to calculate coverage). sample C."
    resources: time_min=180, mem_mb=8000, cpus=4
    benchmark:
        "benchmarks/hisat2_index/tab.tsv"
    conda:
        "envs/hisat2.yaml"
    shell:
        """
        hisat2-build -p 4 {input} {input}
        """


rule mapping_B:
    input:
        reads1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/B_1_dedup_host_removed.fastq.gz",
        reads2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/B_2_dedup_host_removed.fastq.gz",
        input1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.1.ht2",
        input2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.2.ht2",
        input3="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.3.ht2",
        input4="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.4.ht2",
        input5="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.5.ht2",
        input6="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.6.ht2",
        input7="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.7.ht2",
        input8="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa.8.ht2"

    message: "mapping reads on assembly B from files {input.reads1} and {input.reads2}."
    resources: time_min=240, mem_mb=8000, cpus=16
    benchmark:
        "benchmarks/hisat2_align_to_assembly/B.tsv"
    params: 
        db="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa",
        base_out="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/B_1_dedup_host_removed_on_assembly"
    output:
        temp("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/B_1_dedup_host_removed_on_assembly.sorted.bam")
    conda:
        "envs/hisat2.yaml"
    log: "logs/hisat2_mapping_on_assembly_B.log"
    shell:# write only paired-end reads that fail to align concordantly
        """
        hisat2 --summary-file {input.reads1}.hisat_to_assembly.log -x {params.db} -1 {input.reads1} -2 {input.reads2} -k 1 -p 16 --no-discordant -S {params.base_out}.sam 2> {log}
        samtools view -S -b {params.base_out}.sam > {params.base_out}.bam
        samtools sort {params.base_out}.bam -o {output}
        rm {params.base_out}.sam {params.base_out}.bam
        """



rule mapping_C:
    input:
        reads1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_1_dedup_host_removed.fastq.gz",
        reads2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_2_dedup_host_removed.fastq.gz",
        input1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.1.ht2",
        input2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.2.ht2",
        input3="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.3.ht2",
        input4="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.4.ht2",
        input5="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.5.ht2",
        input6="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.6.ht2",
        input7="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.7.ht2",
        input8="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa.8.ht2"

    message: "mapping reads on assembly (C_subset) from files {input.reads1} and {input.reads2}."
    resources: time_min=240, mem_mb=8000, cpus=16
    benchmark:
        "benchmarks/hisat2_align_to_assembly/C_subset.tsv"
    params: 
        db="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa",
        base_out="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_1_dedup_host_removed_on_assembly"
    output:
        temp("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_1_dedup_host_removed_on_assembly.sorted.bam")
    conda:
        "envs/hisat2.yaml"
    log: "logs/hisat2_mapping_on_assembly_C.log"
    shell:
        """
        hisat2 --summary-file {input.reads1}.hisat_to_assembly.log -x {params.db} -1 {input.reads1} -2 {input.reads2} -k 1 -p 16 --no-discordant -S {params.base_out}.sam 2> {log}
        samtools view -S -b {params.base_out}.sam > {params.base_out}.bam 
        samtools sort {params.base_out}.bam -o {output} 
        rm {params.base_out}.sam {params.base_out}.bam 
        """

rule coverage_B:
    input:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/B_1_dedup_host_removed_on_assembly.sorted.bam"
    output:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/B_1_dedup_host_removed_on_assembly.sorted_COV_TABLE"
    resources: time_min=240, mem_mb=8000, cpus=1
    conda:
        "envs/hisat2.yaml"
    benchmark:
        "benchmarks/coverage/B.tsv"
    shell:
        """
        samtools coverage {input} > {output}
        """


rule coverage_C:
    input:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_1_dedup_host_removed_on_assembly.sorted.bam"
    output:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_1_dedup_host_removed_on_assembly.sorted_COV_TABLE"
    resources: time_min=240, mem_mb=8000, cpus=1
    conda:
        "envs/hisat2.yaml"
    benchmark:
        "benchmarks/coverage/C.tsv"
    shell:
        """
        samtools coverage {input} > {output}
        """

from snakemake.utils import R

rule subset_contigs_B:
    message: "Filter out small contigs, lowa quality contigs and contigs with low coverage. Splits the output for blast"
    input:  
        tab="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/B_1_dedup_host_removed_on_assembly.sorted_COV_TABLE",
        contigs="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/final.contigs.fa"
    output: 
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/contigs_to_keep_assembly_B.fa",
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/contigs_filtered_out_assembly_B.fa",
        #"/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/split/file_{sample}.fa"
    conda:
        "envs/biostrings.yaml"
    params: 
        split_dir="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/split"
    script:
        "Rscript_subset_contigs.R"

rule subset_contigs_C:
    message: "Filter out small contigs, low quality contigs and contigs with low coverage. Splits the output for blast"
    input:
        tab="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/trimmed/C_1_dedup_host_removed_on_assembly.sorted_COV_TABLE",
        contigs="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/final.contigs.fa"
    output:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/contigs_to_keep_assembly_C.fa",
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/contigs_filtered_out_assembly_C.fa",
        #"/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/split/file_{sample}.fa"
    conda:
        "envs/biostrings.yaml"
    params: 
        split_dir="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/split"
    script:
        "Rscript_subset_contigs.R"

rule wga_wta_contigs:
    output:
        WGA_contigs="/beegfs/data/varaldi/MACROGEN/WGA-WTA_scaffolds/scaffolds_WGA2.fa",
        WTA_contigs="/beegfs/data/varaldi/MACROGEN/WGA-WTA_scaffolds/scaffolds_WTA2.fa"


rule merge_assemblies:
    message: "merge the assemblies produced from samples B and C and WGA and WTA (independant) assemblies"
    input:
        B_sample="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/contigs_to_keep_assembly_B.fa",
        C_sample="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/contigs_to_keep_assembly_C.fa",
        WGA_contigs="/beegfs/data/varaldi/MACROGEN/WGA-WTA_scaffolds/scaffolds_WGA2.fa",
        WTA_contigs="/beegfs/data/varaldi/MACROGEN/WGA-WTA_scaffolds/scaffolds_WTA2.fa",
    output:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.fasta"
    resources: time_min=120, mem_mb=8000, cpus=8
    conda:
        "envs/flye.yaml"
    log: 
        "logs/mmseqs/flye.log"
    params: 
        out="/beegfs/data/varaldi/MACROGEN/scripts/logs_slurm/flye.out",
        err="/beegfs/data/varaldi/MACROGEN/scripts/logs_slurm/flye.err",
        name="merge"
    benchmark:
        "benchmarks/flye/table.tsv"
    shell:
        """
        hostname
        # modify sequence names to tag their origin
        # sed '/^>/s/\([^ ]*\)\(.*\)/\1_Lh\2/' /beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/contigs_to_keep_assembly_C.fa -i 
        # then we need to fool snakemake, since we modified the file. This can be solved by modifiyig the date of creation of the modified file:
        #touch --date='2021-05-11 17:00:00' /beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/C_subsampling/contigs_to_keep_assembly_C.fa # anti-dating...
        flye --subassemblies {input.B_sample} {input.C_sample} {input.WGA_contigs} {input.WTA_contigs} --out-dir /beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling/ --threads {resources.cpus} --iterations 0 
        echo 'ok merging finished'
        """

rule split_assembly_before_blastx:
    message: "split assembly in files with 300 sequences. Note that braces not used for variable access have to be escaped by repeating them, i.e. {{print $1}}."
    input: 
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.fasta"
    output: 
        expand("/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/split/assembly.fasta-split-{number}.fa", number=n_split)
    params:
        base="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/"
    shell:
        """
        cd {params.base}
        mkdir {params.base}'split'
        awk ' $0~">" {{a++}}{{ b=int(a/300); print $0 > FILENAME"-split-"b".fa"}}' assembly.fasta
        mv assembly.fasta-split-*.fa split
        """


rule blastx_all_contigs_split:
    message: "blastx(mmseqs2) contigs obtained from merging asemblies of sample B (pool samples DNA +cDNAs), sample C( L. heterotoma) and WGA and WTA contigs obtained previously against virus refseq small db."
    input:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/split/assembly.fasta-split-{sample}.fa"
    output:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/blastx/assembly.fasta-split-{sample}.blastx"
    resources: time_min=180, mem_mb=8000, cpus=4
    conda:
        "envs/mmseqs2.yaml"
    log: 
        "logs/mmseqs/full_split-{sample}.log"
    params: 
        db="/beegfs/data/varaldi/MACROGEN/VIRUS_DB/files/refseq_small",
        tmp="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/tmp_mmseqs2-{sample}",
        name="bl-{sample}",
        out="/beegfs/data/varaldi/MACROGEN/scripts/logs_slurm/mmseqs2-split-{sample}.out",
        err="/beegfs/data/varaldi/MACROGEN/scripts/logs_slurm/mmseqs2-split-{sample}.err"
    benchmark:
        "benchmarks/mmseqs/full-split-{sample}.tsv"
    shell:
        """
        hostname
        mmseqs easy-search {input} {params.db}  {output} {params.tmp} --threads 4 --max-seqs 10 -e 0.001
        # clean temp file
        rm -R {params.tmp} 
        echo 'ok blastx finished'
        """



rule hisat2_index_merge:
    input:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.fasta"
    output:
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.1.ht2",
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.2.ht2",
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.3.ht2",
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.4.ht2",
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.5.ht2",
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.6.ht2",
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.7.ht2",
        "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.8.ht2"
    resources: time_min=20, mem_mb=8000, cpus=1
    conda:
        "envs/hisat2.yaml"
    params: 
        base="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly",
        name="hisat2-index",
        out="/beegfs/data/varaldi/MACROGEN/scripts/logs_slurm/hisat2_index_merge.out",
        err="/beegfs/data/varaldi/MACROGEN/scripts/logs_slurm/hisat2_index_merge.err"
    log: 
        "logs/hisat2/create_index_merged.log"
    benchmark:
        "benchmarks/hisat2/create_index_merged.tsv"
    shell:
        """
        hostname 
        hisat2-build -p 1 {input} {params.base}
        echo 'ok index finished'
        """

rule mapping_wgawta_reads:
    message: "map the reads obtained from each sample (after wga or wta) on the db"
    input:
        READS1="/beegfs/data/varaldi/VIROMICS/data/pairs/rename2/trimmed-{number}_R1.fastq",
        READS2="/beegfs/data/varaldi/VIROMICS/data/pairs/rename2/trimmed-{number}_R2.fastq",
        SINGLES="/beegfs/data/varaldi/VIROMICS/data/pairs/rename2/trimmed-{number}_singles.fastq",
        DB="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.1.ht2"

    output: "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_wga_wta_on_merged_assembly/sample_{number}.table.txt"
    resources: time_min=10, mem_mb=750, cpus=1
    conda:
        "envs/hisat2.yaml"
    log: "/beegfs/data/varaldi/MACROGEN/scripts/logs/hisat2-wgawta/{number}.log"
    params: 
        name="map{number}",
        db="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly",
        out="/beegfs/data/varaldi/MACROGEN/scripts/logs_slurm/mapping_wgwta_on_merged_assembly-{number}.out",
        err="/beegfs/data/varaldi/MACROGEN/scripts/logs_slurm/mapping_wgwta_on_merged_assembly-{number}.err",
        out_dir="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_wga_wta_on_merged_assembly",
        base_name="sample_{number}"
    benchmark:
        "benchmarks/hisat2/mapping_wgwta_on_merged_assembly-{number}.tsv"
    shell:
        """
        hostname
        hisat2 --version 2>{log}
        samtools --version 2>{log}
        hisat2 -p 1 --summary-file {params.out_dir}/{params.base_name}.txt -x {params.db} -1 {input.READS1} -2 {input.READS2} -U {input.SINGLES} -S {params.out_dir}/{params.base_name}.sam 
        # convert sam to bam
        samtools view -S -b {params.out_dir}/{params.base_name}.sam > {params.out_dir}/{params.base_name}.bam 
        # sort
        samtools sort {params.out_dir}/{params.base_name}.bam -o {params.out_dir}/{params.base_name}.sorted.bam
        # generate a table of coverage
        samtools coverage {params.out_dir}/{params.base_name}.sorted.bam > {output}
        # clean
        rm {params.out_dir}/{params.base_name}.sam {params.out_dir}/{params.base_name}.bam
        """



rule mapping_reads_on_final_assembly:
    message: "map the reads obtained from each sample (pool or Lh) on the final db"
    input:
        READS1="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_1_dedup.fastq.gz",
        READS2="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/01.RawData/{sample}_2_dedup.fastq.gz",
        DB="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly.1.ht2"

    output: "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_reads_on_merged_assembly/sample_{sample}.table.txt","/beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_reads_on_merged_assembly/sample_{sample}.sorted.bam"
    resources: time_min=240, mem_mb=6000, cpus=16
    conda:
        "envs/hisat2.yaml"
    log: "/beegfs/data/varaldi/MACROGEN/scripts/logs/hisat2-all-reads/{sample}.log"
    params: 
        name="map{sample}",
        db="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B+C_subsampling+WGTA/assembly",
        out="/beegfs/data/varaldi/MACROGEN/scripts/logs_slurm/mapping_reads_on_merged_assembly-{sample}.out",
        err="/beegfs/data/varaldi/MACROGEN/scripts/logs_slurm/mapping_reads_on_merged_assembly-{sample}.err",
        out_dir="/beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_reads_on_merged_assembly",
        base_name="sample_{sample}"
    benchmark:
        "benchmarks/hisat2/mapping_reads_on_merged_assembly-{sample}.tsv"
    shell:
        """
        hostname
        hisat2 -p 16 --summary-file {params.out_dir}/{params.base_name}.txt -x {params.db} -1 {input.READS1} -2 {input.READS2} -S {params.out_dir}/{params.base_name}.sam 
        # convert sam to bam
        samtools view -S -b {params.out_dir}/{params.base_name}.sam > {params.out_dir}/{params.base_name}.bam 
        # sort
        samtools sort {params.out_dir}/{params.base_name}.bam -o {params.out_dir}/{params.base_name}.sorted.bam
        # generate a table of coverage
        samtools coverage {params.out_dir}/{params.base_name}.sorted.bam > {output}
        # clean
        rm {params.out_dir}/{params.base_name}.sam {params.out_dir}/{params.base_name}.bam
        """


rule coverage_depth:
    message: "calculate coverage depth for contigs of interest (FINAL_contigs)"
    input:"/beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_reads_on_merged_assembly/sample_B.sorted.bam","/beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_reads_on_merged_assembly/contigs_of_interest_wga_wta.txt"
    output: "/beegfs/data/varaldi/MACROGEN/2102KNO-0022/mapping_reads_on_merged_assembly/contigs_of_interest_wga_wta.depth.txt"
    shell:
        """
        sh coverage_depth.sh
        """


