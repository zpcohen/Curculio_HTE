### this script was utilized to generate reliable metrics regarding size and position of nodes per species. It was run on each chromosome graph. 


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

dir="/mnt/archive/groups/degiorgio_group/zpc/Curculio_catChromos/LowerPI_Cglan5_pggb"
graph=$(ls "$dir"/*.og | head -n1)

RUN="/mnt/beegfs/home/zcohen5/micromamba/bin/micromamba run -p /mnt/archive/groups/degiorgio_group/zpc/.mamba/envs/pggb-074"
ODGI="/mnt/archive/groups/degiorgio_group/zpc/.mamba/envs/pggb-074/bin/odgi"

# All positions across all paths
$RUN odgi position -i "$graph" --all-positions > "$dir/LowerPI_Cglan5_pggb.tsv"

# Node lengths (don’t append; overwrite)
$RUN odgi view -i "$graph" -g | awk -F'\t' '$1=="S"{
  len="";
  if ($3!="*") len=length($3);            # sequence present
  for(i=4;i<=NF;i++) if($i~/^LN:i:/){     # or LN:i:<len> tag
    split($i,a,":"); len=a[3]
  }
  print $2, len
}' > "$dir/LowerPI_Cglan5_pggb_lengths.tsv"

# Keep nodes >=100 bp
awk '$2 >= 100' "$dir/LowerPI_Cglan5_pggb_lengths.tsv" > "$dir/100bp_LowerPI_Cglan5_pggb_node_lengths.tsv"

# Join node lengths onto positions (fix filename + consistency)
LC_ALL=C awk -F'[[:space:]]+' -v OFS='\t' '
  NR==FNR { len[$1]=$2; next }           # lengths already filtered >=100
  FNR==1   { print $0, "node_len_bp"; next }
  ($2 in len) { print $0, len[$2] }
' "$dir/100bp_LowerPI_Cglan5_pggb_node_lengths.tsv" "$dir/LowerPI_Cglan5_pggb.tsv" > "$dir/nodes_LowerPI_Cglan5_pggb_100bp.tsv"

# Build comma-separated list of reference paths starting with "Cg" (fix assignment + stray quote)
REF=$(grep -E '^Cg' "$dir/nodes_LowerPI_Cglan5_pggb_100bp.tsv" | cut -f1 | sort -u | paste -sd, -)

# Positions constrained to reference paths
$RUN odgi position -i "$graph" --all-positions --ref-path "$REF" > "$dir/LowerPI_Cglan5_pggb_CgREFPATH.tsv"

# Node ID list (skip header)
awk 'NR>1{print $2}' "$dir/nodes_LowerPI_Cglan5_pggb_100bp.tsv" > "$dir/nodes_LowerPI_Cglan5_pggb_100bp.list"

# Project node positions onto the reference (don’t append; overwrite)
$ODGI position -i "$graph" -G "$dir/nodes_LowerPI_Cglan5_pggb_100bp.list" -r "$REF" -d 15000000 -w 15000000 > "$dir/LowerPI_Cglan5_pggb_onREF.tsv"

### GENERATE R FILES (simplify the comma split)
awk -F',' '{print $1}' "$dir/LowerPI_Cglan5_pggb_onREF.tsv" > "$dir/nodeID.col2"
awk -F',' '{print $2}' "$dir/LowerPI_Cglan5_pggb_onREF.tsv" > "$dir/Ref_start.col2"
awk '{print $1}'  "$dir/nodes_LowerPI_Cglan5_pggb_100bp.tsv" > "$dir/Source_col1.list"
awk '{print $4}'  "$dir/nodes_LowerPI_Cglan5_pggb_100bp.tsv" > "$dir/node_length.col4"
awk -F',' '{print $1}' "$dir/LowerPI_Cglan5_pggb_onREF.tsv" > "$dir/target.path.pos_col5"

paste -d$'\t' "$dir/Source_col1.list" "$dir/nodeID.col2" "$dir/Ref_start.col2" "$dir/node_length.col4" "$dir/target.path.pos_col5" > "$dir/Final_LowerPI_Cglan5_pggb_SVs.tsv"
