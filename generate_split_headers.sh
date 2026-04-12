grep '^>' Split_FastasALL.fa | sed 's/^>//' | awk -v OFS='\t' -v w=5000 '
{
    n = split($0, a, "_")
    parent = a[1]
    for (i=2; i<=n-2; i++) parent = parent "_" a[i]
    start = a[n-1]
    end   = a[n]
    s1 = start - 1
    e1 = start - 1 + w
    if (s1 < 0) s1 = 0
    print parent, s1, e1, "chunk_start"
    s2 = end - w
    if (s2 < 0) s2 = 0
    e2 = end
    print parent, s2, e2, "chunk_end"
}' >> Split_Fastas_Headers.bed
