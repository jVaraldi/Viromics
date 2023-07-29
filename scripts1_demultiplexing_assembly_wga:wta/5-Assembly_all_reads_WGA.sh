#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/beegfs/data/varaldi/VIROMICS/data/assemblies/ALL_together/assemblyWGA.log.out
#SBATCH --error=/beegfs/data/varaldi/VIROMICS/data/assemblies/ALL_together/assemblyWGA.log.error

cd /beegfs/data/varaldi/VIROMICS/data/assemblies/ALL_together

# assembling all reads from the metagenomic survey (WGA libs only)

readsF=/beegfs/data/varaldi/VIROMICS/data/pairs/WGA_ALL_R1.fastq
readsR=/beegfs/data/varaldi/VIROMICS/data/pairs/WGA_ALL_R2.fastq
singles=/beegfs/data/varaldi/VIROMICS/data/pairs/WGA_ALL_singles.fastq

#~/TOOLS/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 $readsF -2 $readsR -r
$singles -o out_wga
