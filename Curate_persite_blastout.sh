find BLAST_all_Ccar*2/*out_new*/ -type f -name "*blast.out" -exec awk '
BEGIN{FS=OFS="\t"}
$4 >= 1000 && $11 <= 1e-10 && $3 >= 70 {
  s=($9<$10?$9:$10); e=($9<$10?$10:$9)
  # chr start end qid evalue bitscore length
  print $2, s, e, $1, $11, $12, $4
}' {} + \
| sort -k1,1 -k2,2n \
| awk '
BEGIN{FS=OFS="\t"}

function better(ev, bs, bev, bbs) {
  if (ev < bev) return 1
  if (ev == bev && bs > bbs) return 1
  return 0
}

NR==1 {
  chr=$1; rs=$2; re=$3
  n=1
  best_id=$4; best_ev=$5+0; best_bs=$6+0; best_len=$7+0
  next
}

{
  c=$1; s=$2; e=$3; id=$4; ev=$5+0; bs=$6+0; len=$7+0

  if (c==chr && s<=re) {
    # overlap: extend locus
    if (e>re) re=e
    n++
    if (better(ev,bs,best_ev,best_bs)) { best_id=id; best_ev=ev; best_bs=bs; best_len=len }
  } else {
    print chr, rs, re, n, best_id, best_ev, best_bs, best_len
    chr=c; rs=s; re=e
    n=1
    best_id=id; best_ev=ev; best_bs=bs; best_len=len
  }
}

END {
  if (NR) print chr, rs, re, n, best_id, best_ev, best_bs, best_len
}'
