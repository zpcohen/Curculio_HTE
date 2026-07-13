#!/usr/bin/env python3
"""
blast_divergence.py

Convert BLAST tabular (outfmt 6) percent identity into evolutionary distances.

Pipeline per alignment unit:
    pident  ->  identities / aligned columns
            ->  p-distance (gapped or gap-excluded basis)
            ->  substitution-model-corrected distance

Nucleotide models: p-distance, Jukes-Cantor (1969).
Protein models:    p-distance, Poisson, Kimura (1983) empirical.

Because NCBI BLAST computes pident over an alignment length that includes gap
columns, two p-distances are always reported:
    p_gapped : 1 - pident/100, so indel columns count as differences.
    p_nogap  : mismatches / (identities + mismatches), gap columns removed.
Substitution models model substitutions and not indels, so the default model
input is p_nogap. If you want the BLAST-native behaviour, then pass
--basis gapped.

Multiple HSPs per query-subject pair are handled by --aggregate:
    none : one row per HSP (default, fully transparent).
    best : single best-bitscore HSP per query-subject.
    pool : length-weighted pooling, but only within a subject locus. HSPs that
           map to disjoint subject regions (separated by more than --locus-gap,
           or on opposite strands) are treated as separate loci and pooled
           separately, so a duplicated gene does not get averaged across copies.

Standard library only. No pandas, no numpy.
"""

import argparse
import csv
import math
import sys
from dataclasses import dataclass

NAN = float("nan")

DEFAULT_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]

# --------------------------------------------------------------------------- #
# Distance models. Each maps a p-distance to a corrected distance and returns
# NaN when p is at or beyond the model ceiling (saturation).
# --------------------------------------------------------------------------- #

def d_pdistance(p):
    return p

def d_jc69(p):
    """Jukes-Cantor 1969 nucleotide distance. Ceiling at p = 0.75."""
    arg = 1.0 - (4.0 / 3.0) * p
    if arg <= 0.0:
        return NAN
    return -0.75 * math.log(arg)

def d_poisson(p):
    """Poisson-corrected protein distance. Ceiling at p = 1."""
    arg = 1.0 - p
    if arg <= 0.0:
        return NAN
    return -math.log(arg)

def d_kimura_protein(p):
    """Kimura 1983 empirical protein distance. Ceiling near p = 0.8541."""
    arg = 1.0 - p - 0.2 * p * p
    if arg <= 0.0:
        return NAN
    return -math.log(arg)

NUCL_MODELS = {"p_distance": d_pdistance, "jc69": d_jc69}
PROT_MODELS = {"p_distance": d_pdistance, "poisson": d_poisson, "kimura": d_kimura_protein}

# --------------------------------------------------------------------------- #
# Records
# --------------------------------------------------------------------------- #

@dataclass
class Hsp:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    sstart: int
    send: int
    bitscore: float

    def __post_init__(self):
        # Derive identity count from pident and length, then recover total gap
        # characters as the residual. This sidesteps the fact that outfmt 6
        # reports gap openings and not total gap length.
        self.identities = int(round(self.pident / 100.0 * self.length))
        residual = self.length - self.identities - self.mismatch
        self.gaps = residual if residual > 0 else 0

    @property
    def slo(self):
        return min(self.sstart, self.send)

    @property
    def shi(self):
        return max(self.sstart, self.send)

    @property
    def strand(self):
        return "+" if self.send >= self.sstart else "-"


@dataclass
class AlnUnit:
    """One row of output: a single HSP, a best HSP, or a pooled locus."""
    qseqid: str
    sseqid: str
    n_hsp: int
    aln_length: int
    identities: int
    mismatch: int
    gaps: int
    slo: int
    shi: int
    strand: str
    bitscore: float

    @classmethod
    def from_hsps(cls, hsps):
        aln_length = sum(h.length for h in hsps)
        identities = sum(h.identities for h in hsps)
        mismatch = sum(h.mismatch for h in hsps)
        gaps = sum(h.gaps for h in hsps)
        strands = {h.strand for h in hsps}
        strand = next(iter(strands)) if len(strands) == 1 else "mixed"
        return cls(
            qseqid=hsps[0].qseqid,
            sseqid=hsps[0].sseqid,
            n_hsp=len(hsps),
            aln_length=aln_length,
            identities=identities,
            mismatch=mismatch,
            gaps=gaps,
            slo=min(h.slo for h in hsps),
            shi=max(h.shi for h in hsps),
            strand=strand,
            bitscore=sum(h.bitscore for h in hsps),
        )

    @property
    def pident(self):
        return 100.0 * self.identities / self.aln_length if self.aln_length else NAN

    @property
    def p_gapped(self):
        return 1.0 - self.identities / self.aln_length if self.aln_length else NAN

    @property
    def p_nogap(self):
        ungapped = self.identities + self.mismatch
        return self.mismatch / ungapped if ungapped > 0 else NAN

    def p_on(self, basis):
        return self.p_nogap if basis == "nogap" else self.p_gapped

# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

def resolve_indices(columns):
    """Map required field names to positions in a user-declared column list."""
    required = ["qseqid", "sseqid", "pident", "length", "mismatch",
                "gapopen", "sstart", "send", "bitscore"]
    idx = {}
    for name in required:
        if name not in columns:
            sys.exit(f"error: required column '{name}' not present in --columns")
        idx[name] = columns.index(name)
    return idx


def read_hsps(path, columns, min_aln_length):
    idx = resolve_indices(columns)
    handle = sys.stdin if path == "-" else open(path, newline="")
    hsps = []
    try:
        reader = csv.reader(handle, delimiter="\t")
        for lineno, row in enumerate(reader, 1):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < len(columns):
                sys.exit(f"error: line {lineno} has {len(row)} fields, "
                         f"expected {len(columns)}")
            try:
                h = Hsp(
                    qseqid=row[idx["qseqid"]],
                    sseqid=row[idx["sseqid"]],
                    pident=float(row[idx["pident"]]),
                    length=int(row[idx["length"]]),
                    mismatch=int(row[idx["mismatch"]]),
                    gapopen=int(row[idx["gapopen"]]),
                    sstart=int(row[idx["sstart"]]),
                    send=int(row[idx["send"]]),
                    bitscore=float(row[idx["bitscore"]]),
                )
            except ValueError as exc:
                sys.exit(f"error: could not parse line {lineno}: {exc}")
            if h.length >= min_aln_length:
                hsps.append(h)
    finally:
        if handle is not sys.stdin:
            handle.close()
    return hsps

# --------------------------------------------------------------------------- #
# Aggregation
# --------------------------------------------------------------------------- #

def group_by_pair(hsps):
    groups = {}
    for h in hsps:
        groups.setdefault((h.qseqid, h.sseqid), []).append(h)
    return groups


def cluster_loci(hsps, locus_gap):
    """Partition a set of same-pair HSPs into subject loci.

    HSPs are merged into one locus while they share a strand and their subject
    intervals are separated by no more than locus_gap. A change of strand or a
    larger gap starts a new locus. This keeps split alignments of one gene
    together while separating duplicated copies.
    """
    loci = []
    for strand in ("+", "-"):
        same = sorted((h for h in hsps if h.strand == strand), key=lambda h: h.slo)
        current = []
        current_hi = None
        for h in same:
            if current and h.slo <= current_hi + locus_gap:
                current.append(h)
                current_hi = max(current_hi, h.shi)
            else:
                if current:
                    loci.append(current)
                current = [h]
                current_hi = h.shi
        if current:
            loci.append(current)
    return loci


def build_units(hsps, aggregate, locus_gap):
    if aggregate == "none":
        return [AlnUnit.from_hsps([h]) for h in hsps]

    units = []
    for pair_hsps in group_by_pair(hsps).values():
        if aggregate == "best":
            best = max(pair_hsps, key=lambda h: h.bitscore)
            units.append(AlnUnit.from_hsps([best]))
        elif aggregate == "pool":
            for locus in cluster_loci(pair_hsps, locus_gap):
                units.append(AlnUnit.from_hsps(locus))
        else:
            sys.exit(f"error: unknown aggregate mode '{aggregate}'")
    return units

# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #

def fmt(x, places=6):
    if x is None:
        return "NA"
    if isinstance(x, float) and math.isnan(x):
        return "NA"
    if isinstance(x, float):
        return f"{x:.{places}f}"
    return str(x)


def write_table(units, models, basis, out_handle):
    model_names = list(models.keys())
    header = (
        ["qseqid", "sseqid", "n_hsp", "aln_length", "identities", "mismatch",
         "gaps", "pident", "p_gapped", "p_nogap", "basis"]
        + [f"d_{name}" for name in model_names]
        + ["s_start", "s_end", "strand", "bitscore", "flag"]
    )
    writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
    writer.writerow(header)

    for u in sorted(units, key=lambda u: (u.qseqid, u.slo)):
        p = u.p_on(basis)
        dvals = {name: fn(p) for name, fn in models.items()}
        saturated = any(isinstance(v, float) and math.isnan(v)
                        for name, v in dvals.items() if name != "p_distance")
        flag = "saturated" if saturated else ""
        row = (
            [u.qseqid, u.sseqid, u.n_hsp, u.aln_length, u.identities,
             u.mismatch, u.gaps, fmt(u.pident, 3),
             fmt(u.p_gapped), fmt(u.p_nogap), basis]
            + [fmt(dvals[name]) for name in model_names]
            + [u.slo, u.shi, u.strand, fmt(u.bitscore, 1), flag]
        )
        writer.writerow(row)

# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Convert BLAST outfmt 6 percent identity to evolutionary distances.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("blast_table", help="BLAST outfmt 6 file, or - for stdin")
    p.add_argument("-o", "--out", default="-",
                   help="output TSV path, or - for stdout")
    p.add_argument("--seqtype", choices=["nucl", "prot"], default="nucl",
                   help="nucleotide (JC69) or protein (Poisson, Kimura) models")
    p.add_argument("--basis", choices=["nogap", "gapped"], default="nogap",
                   help="p-distance fed to substitution models")
    p.add_argument("--aggregate", choices=["none", "best", "pool"], default="none",
                   help="per-HSP, best HSP per pair, or pooled per subject locus")
    p.add_argument("--locus-gap", type=int, default=5000,
                   help="max subject-coordinate gap (bp) merged into one locus "
                        "when --aggregate pool")
    p.add_argument("--min-aln-length", type=int, default=0,
                   help="drop HSPs shorter than this many columns")
    p.add_argument("--columns", default=",".join(DEFAULT_COLUMNS),
                   help="comma-separated outfmt 6 column names in file order")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    columns = [c.strip() for c in args.columns.split(",")]
    models = NUCL_MODELS if args.seqtype == "nucl" else PROT_MODELS

    hsps = read_hsps(args.blast_table, columns, args.min_aln_length)
    if not hsps:
        sys.exit("error: no HSPs parsed")

    units = build_units(hsps, args.aggregate, args.locus_gap)

    out_handle = sys.stdout if args.out == "-" else open(args.out, "w", newline="")
    try:
        write_table(units, models, args.basis, out_handle)
    finally:
        if out_handle is not sys.stdout:
            out_handle.close()

    sys.stderr.write(
        f"parsed {len(hsps)} HSPs, wrote {len(units)} rows "
        f"[seqtype={args.seqtype}, basis={args.basis}, aggregate={args.aggregate}]\n"
    )


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # Downstream reader (head, less) closed the pipe. Exit quietly.
        try:
            sys.stdout.close()
        except Exception:
            pass
        sys.exit(0)
