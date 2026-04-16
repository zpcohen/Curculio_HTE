#!/bin/bash
#SBATCH -J pggb
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=100G
#SBATCH -p longq7-eng
#SBATCH -t 48:00:00
#SBATCH -o logs/C6_%x-%j.out
set -euo pipefail


RUN="/mnt/beegfs/home/zcohen5/micromamba/bin/micromamba run -p /mnt/archive/groups/degiorgio_group/zpc/.mamba/envs/pggb-074"

$RUN pggb -i Cg6Cc7Cn7_new.fasta -o /mnt/archive/groups/degiorgio_group/zpc/Curculio_catChromos/LowerPI_Cglan6_pggb \
	-n 3 \
 	-t 32 \
	-p 50 \
	-r \
	-s 1000 \
	-V 'Cglandium:1000'
