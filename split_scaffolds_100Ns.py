#!/usr/bin/env python3
import re
import sys

if len(sys.argv) < 3:
    sys.stderr.write("usage: split_fasta_on_N.py input.fa min_N_run > output.fa\n")
    sys.exit(1)

fasta = sys.argv[1]
min_n = int(sys.argv[2])

pat = re.compile(r"[Nn]{%d,}" % min_n)

def write_chunks(name, seq):
    seq = "".join(seq)
    parts = pat.split(seq)
    coords = []
    last = 0
    for m in pat.finditer(seq):
        coords.append((last, m.start()))
        last = m.end()
    coords.append((last, len(seq)))

    chunk_id = 0
    for part, (start, end) in zip(parts, coords):
        if not part:
            continue
        chunk_id += 1
        header = f">{name}_chunk{chunk_id}_{start+1}_{end}"
        print(header)
        for i in range(0, len(part), 80):
            print(part[i:i+80])

name = None
seq = []

with open(fasta) as fh:
    for line in fh:
        line = line.rstrip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                write_chunks(name, seq)
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)

if name is not None:
    write_chunks(name, seq)

### then run generate_split_headers.sh
