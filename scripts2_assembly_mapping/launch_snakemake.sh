#!/bin/bash
conda activate snakemake
snakemake --profile slurm --rerun-incomplete -p  /beegfs/data/varaldi/MACROGEN/2102KNO-0022/assembly/B/blastx/file_{1,2,3,4,5}.blastx
