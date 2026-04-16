#!/bin/bash
#SBATCH -J pggb
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=100G
#SBATCH -p longq7-eng
#SBATCH -t 48:00:00
#SBATCH -o logs/%x-%j.out
set -euo pipefail


#!/bin/bash
set -euo pipefail


RUN="/mnt/beegfs/home/zcohen5/micromamba/bin/micromamba run -p /mnt/archive/groups/degiorgio_group/zpc/.mamba/envs/pggb-074"

$RUN odgi position -i CglanCcarCnanu3_new.fasta.2afcd0e.11fba48.33b105f.smooth.final.og -G Ccaryae#3#Chrom3_list -r Ccaryae#3#Chrom3 > Ccaryae_C3_90pc_node_start2.tsv
