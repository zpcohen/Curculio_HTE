#!/bin/bash
#SBATCH -J earlgrey_RUN_Ccaryae5
#SBATCH -p longq7-eng
##SBATCH --cpus-per-task=72
#SBATCH -N 1
#SBATCH -n 72
#SBATCH --mem=800G
#SBATCH -t 168:00:00
#SBATCH -o logs/RUN_Cc2Feb20%x-%j.out
#SBATCH -e logs/RUN_Cc2Feb20%x-%j.err

set -euo pipefail

export MAMBA_ROOT_PREFIX=/mnt/archive/groups/degiorgio_group/zpc/.mamba

$HOME/micromamba/bin/micromamba run -p "$MAMBA_ROOT_PREFIX/envs/earlgrey" earlGrey -t 72 -g /mnt/archive/groups/degiorgio_group/zpc/Curculio_catChromos/Ccaryae_Feb24.fasta -s Ccaryae-UM -o earlgrey_TE_output/Cc_C2_out_Feb24/
