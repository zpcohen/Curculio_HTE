awk '
  BEGIN { OFS="\t" }

  # PASS 1: read odgi position output and build node -> (chrom,start)
  FNR==NR {
    if ($0 ~ /^#|^$/) next

    # split line into tab fields
    n = split($0, f, "\t")

    # f[1] = "1597564,0,+"
    # f[2] = "Cnanu#1#Chrom1,28850,+"
    split(f[1], a, ","); node = a[1]
    split(f[2], b, ","); chrom = b[1]; start = b[2]

    chr[node] = chrom
    pos[node] = start
    next
  }

  # PASS 2: read nodes.tsv and emit intervals
  {
    size = $4
    node = $5

    if (!(node in pos)) next

    s = pos[node] + 0
    e = s + size

    print chr[node], s, e
  }
' Cnanu_C3_90pc_node_start2.tsv ../CnanuMar6_pangenomenodes_90pc.tsv 
