#!/usr/bin/env python3
"""
blast_divergence_weighted.py

Convert BLAST tabular output (outfmt 6) into nucleotide or protein
evolutionary distances.

For nucleotide data, the script reports:
    - p-distance
    - Jukes-Cantor 1969 (JC69) distance

For protein data, the script reports:
    - p-distance
    - Poisson-corrected distance
    - Kimura (1983) empirical distance

Two alignment bases are supported:
    nogap  : mismatches / (identities + mismatches)
    gapped : 1 - identities / alignment_length

Substitution models describe substitutions rather than indels, so ``nogap`` is
the default.

Multiple HSPs per query-subject pair can be handled using ``--aggregate``:

    none
        One output row per HSP.

    best
        The highest-bitscore HSP for each query-subject pair.

    pool
        HSPs are clustered into subject loci and pooled within each locus.
        Opposite strands or HSPs separated by more than ``--locus-gap`` are
        treated as separate loci.

For pooled loci, two corrected-distance summaries are reported:

    d_<model>
        The model correction applied to the pooled site-level p-distance.
        This is the preferred length-weighted regional estimate because every
        aligned site contributes once to the pooled mismatch proportion.

    weighted_mean_d_<model>_hsp
        The arithmetic mean of the per-HSP corrected distances, weighted by
        the number of sites used for the selected basis. This is retained as a
        diagnostic because averaging corrected distances is not identical to
        correcting the pooled p-distance.

Standard library only; no pandas or NumPy required.
"""

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

NAN = float("nan")

DEFAULT_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


# --------------------------------------------------------------------------- #
# Distance models
# --------------------------------------------------------------------------- #

def d_pdistance(p: float) -> float:
    """Uncorrected p-distance."""
    return p


def d_jc69(p: float) -> float:
    """
    Jukes-Cantor 1969 nucleotide distance.

    The correction is undefined when p >= 0.75 because the logarithm argument
    becomes zero or negative.
    """
    if math.isnan(p) or p < 0.0:
        return NAN

    arg = 1.0 - (4.0 / 3.0) * p
    if arg <= 0.0:
        return NAN

    return -0.75 * math.log(arg)


def d_poisson(p: float) -> float:
    """Poisson-corrected protein distance."""
    if math.isnan(p) or p < 0.0:
        return NAN

    arg = 1.0 - p
    if arg <= 0.0:
        return NAN

    return -math.log(arg)


def d_kimura_protein(p: float) -> float:
    """Kimura 1983 empirical protein distance."""
    if math.isnan(p) or p < 0.0:
        return NAN

    arg = 1.0 - p - 0.2 * p * p
    if arg <= 0.0:
        return NAN

    return -math.log(arg)


NUCL_MODELS: Dict[str, Callable[[float], float]] = {
    "p_distance": d_pdistance,
    "jc69": d_jc69,
}

PROT_MODELS: Dict[str, Callable[[float], float]] = {
    "p_distance": d_pdistance,
    "poisson": d_poisson,
    "kimura": d_kimura_protein,
}


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

    def __post_init__(self) -> None:
        if self.length < 0:
            raise ValueError("alignment length cannot be negative")
        if self.mismatch < 0:
            raise ValueError("mismatch count cannot be negative")
        if not 0.0 <= self.pident <= 100.0:
            raise ValueError("pident must be between 0 and 100")

        # BLAST outfmt 6 reports gap openings, not total gapped columns.
        # Recover identity count from pident and alignment length, then treat
        # the remaining columns as gaps.
        self.identities = int(round(self.pident / 100.0 * self.length))
        residual = self.length - self.identities - self.mismatch
        self.gaps = max(residual, 0)

    @property
    def slo(self) -> int:
        return min(self.sstart, self.send)

    @property
    def shi(self) -> int:
        return max(self.sstart, self.send)

    @property
    def strand(self) -> str:
        return "+" if self.send >= self.sstart else "-"

    @property
    def p_gapped(self) -> float:
        if self.length == 0:
            return NAN
        return 1.0 - self.identities / self.length

    @property
    def p_nogap(self) -> float:
        ungapped = self.identities + self.mismatch
        if ungapped == 0:
            return NAN
        return self.mismatch / ungapped

    def p_on(self, basis: str) -> float:
        return self.p_nogap if basis == "nogap" else self.p_gapped

    def weight_on(self, basis: str) -> int:
        """
        Number of sites contributing to the selected p-distance.

        For nogap distances, use identities + mismatches.
        For gapped distances, use the complete BLAST alignment length.
        """
        if basis == "nogap":
            return self.identities + self.mismatch
        return self.length


@dataclass
class AlnUnit:
    """One output unit: one HSP, a best HSP, or a pooled subject locus."""

    qseqid: str
    sseqid: str
    hsps: Sequence[Hsp]

    @property
    def n_hsp(self) -> int:
        return len(self.hsps)

    @property
    def aln_length(self) -> int:
        return sum(h.length for h in self.hsps)

    @property
    def identities(self) -> int:
        return sum(h.identities for h in self.hsps)

    @property
    def mismatch(self) -> int:
        return sum(h.mismatch for h in self.hsps)

    @property
    def gaps(self) -> int:
        return sum(h.gaps for h in self.hsps)

    @property
    def slo(self) -> int:
        return min(h.slo for h in self.hsps)

    @property
    def shi(self) -> int:
        return max(h.shi for h in self.hsps)

    @property
    def strand(self) -> str:
        strands = {h.strand for h in self.hsps}
        return next(iter(strands)) if len(strands) == 1 else "mixed"

    @property
    def bitscore(self) -> float:
        return sum(h.bitscore for h in self.hsps)

    @property
    def pident(self) -> float:
        if self.aln_length == 0:
            return NAN
        return 100.0 * self.identities / self.aln_length

    @property
    def p_gapped(self) -> float:
        if self.aln_length == 0:
            return NAN
        return 1.0 - self.identities / self.aln_length

    @property
    def p_nogap(self) -> float:
        ungapped = self.identities + self.mismatch
        if ungapped == 0:
            return NAN
        return self.mismatch / ungapped

    def p_on(self, basis: str) -> float:
        """
        Pooled site-level p-distance for the complete alignment unit.

        Applying the substitution correction to this value yields the preferred
        length-weighted corrected distance.
        """
        return self.p_nogap if basis == "nogap" else self.p_gapped

    def weighted_mean_model_distance(
        self,
        model: Callable[[float], float],
        basis: str,
    ) -> float:
        """
        Length-weighted mean of per-HSP model-corrected distances.

        HSPs for which the model is saturated or undefined are omitted. If all
        HSPs are undefined, return NaN.
        """
        numerator = 0.0
        denominator = 0

        for hsp in self.hsps:
            weight = hsp.weight_on(basis)
            if weight <= 0:
                continue

            distance = model(hsp.p_on(basis))
            if math.isnan(distance):
                continue

            numerator += weight * distance
            denominator += weight

        if denominator == 0:
            return NAN

        return numerator / denominator


# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #

def resolve_indices(columns: Sequence[str]) -> Dict[str, int]:
    """Map required BLAST field names to positions."""
    required = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "sstart",
        "send",
        "bitscore",
    ]

    missing = [name for name in required if name not in columns]
    if missing:
        sys.exit(
            "error: required columns absent from --columns: "
            + ", ".join(missing)
        )

    return {name: columns.index(name) for name in required}


def read_hsps(
    path: str,
    columns: Sequence[str],
    min_aln_length: int,
) -> List[Hsp]:
    """Read BLAST outfmt 6 records."""
    idx = resolve_indices(columns)
    handle = sys.stdin if path == "-" else open(path, newline="")
    hsps: List[Hsp] = []

    try:
        reader = csv.reader(handle, delimiter="\t")

        for lineno, row in enumerate(reader, 1):
            if not row or row[0].startswith("#"):
                continue

            if len(row) < len(columns):
                sys.exit(
                    f"error: line {lineno} has {len(row)} fields; "
                    f"expected at least {len(columns)}"
                )

            try:
                hsp = Hsp(
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

            if hsp.length >= min_aln_length:
                hsps.append(hsp)

    finally:
        if handle is not sys.stdin:
            handle.close()

    return hsps


# --------------------------------------------------------------------------- #
# Aggregation
# --------------------------------------------------------------------------- #

def group_by_pair(hsps: Iterable[Hsp]) -> Dict[Tuple[str, str], List[Hsp]]:
    groups: Dict[Tuple[str, str], List[Hsp]] = {}

    for hsp in hsps:
        groups.setdefault((hsp.qseqid, hsp.sseqid), []).append(hsp)

    return groups


def cluster_loci(hsps: Sequence[Hsp], locus_gap: int) -> List[List[Hsp]]:
    """
    Partition same-query, same-subject HSPs into subject-coordinate loci.

    HSPs are combined while they:
      1. occur on the same strand; and
      2. overlap or are separated by no more than ``locus_gap`` bases.
    """
    loci: List[List[Hsp]] = []

    for strand in ("+", "-"):
        strand_hsps = sorted(
            (h for h in hsps if h.strand == strand),
            key=lambda h: h.slo,
        )

        current: List[Hsp] = []
        current_hi = None

        for hsp in strand_hsps:
            if (
                current
                and current_hi is not None
                and hsp.slo <= current_hi + locus_gap
            ):
                current.append(hsp)
                current_hi = max(current_hi, hsp.shi)
            else:
                if current:
                    loci.append(current)

                current = [hsp]
                current_hi = hsp.shi

        if current:
            loci.append(current)

    return loci


def make_unit(hsps: Sequence[Hsp]) -> AlnUnit:
    if not hsps:
        raise ValueError("cannot create an alignment unit from zero HSPs")

    return AlnUnit(
        qseqid=hsps[0].qseqid,
        sseqid=hsps[0].sseqid,
        hsps=tuple(hsps),
    )


def build_units(
    hsps: Sequence[Hsp],
    aggregate: str,
    locus_gap: int,
) -> List[AlnUnit]:
    if aggregate == "none":
        return [make_unit([hsp]) for hsp in hsps]

    units: List[AlnUnit] = []

    for pair_hsps in group_by_pair(hsps).values():
        if aggregate == "best":
            best = max(pair_hsps, key=lambda h: h.bitscore)
            units.append(make_unit([best]))

        elif aggregate == "pool":
            for locus_hsps in cluster_loci(pair_hsps, locus_gap):
                units.append(make_unit(locus_hsps))

        else:
            sys.exit(f"error: unknown aggregate mode '{aggregate}'")

    return units


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #

def fmt(value, places: int = 6) -> str:
    if value is None:
        return "NA"

    if isinstance(value, float) and math.isnan(value):
        return "NA"

    if isinstance(value, float):
        return f"{value:.{places}f}"

    return str(value)


def write_table(
    units: Sequence[AlnUnit],
    models: Dict[str, Callable[[float], float]],
    basis: str,
    out_handle,
) -> None:
    model_names = list(models.keys())

    header = [
        "qseqid",
        "sseqid",
        "n_hsp",
        "aln_length",
        "ungapped_length",
        "identities",
        "mismatch",
        "gaps",
        "pident",
        "p_gapped",
        "p_nogap",
        "basis",
    ]

    # d_<model> is the model correction applied to the pooled p-distance.
    header.extend(f"d_{name}" for name in model_names)

    # Diagnostic: weighted mean of independently corrected HSP distances.
    header.extend(
        f"weighted_mean_d_{name}_hsp"
        for name in model_names
    )

    header.extend(
        [
            "s_start",
            "s_end",
            "strand",
            "bitscore",
            "flag",
        ]
    )

    writer = csv.writer(
        out_handle,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writerow(header)

    for unit in sorted(
        units,
        key=lambda u: (u.qseqid, u.sseqid, u.slo, u.shi),
    ):
        pooled_p = unit.p_on(basis)

        pooled_distances = {
            name: model(pooled_p)
            for name, model in models.items()
        }

        weighted_hsp_distances = {
            name: unit.weighted_mean_model_distance(model, basis)
            for name, model in models.items()
        }

        saturated_models = [
            name
            for name, value in pooled_distances.items()
            if name != "p_distance"
            and isinstance(value, float)
            and math.isnan(value)
        ]

        flag = (
            "saturated:" + ",".join(saturated_models)
            if saturated_models
            else ""
        )

        row = [
            unit.qseqid,
            unit.sseqid,
            unit.n_hsp,
            unit.aln_length,
            unit.identities + unit.mismatch,
            unit.identities,
            unit.mismatch,
            unit.gaps,
            fmt(unit.pident, 3),
            fmt(unit.p_gapped),
            fmt(unit.p_nogap),
            basis,
        ]

        row.extend(fmt(pooled_distances[name]) for name in model_names)
        row.extend(
            fmt(weighted_hsp_distances[name])
            for name in model_names
        )

        row.extend(
            [
                unit.slo,
                unit.shi,
                unit.strand,
                fmt(unit.bitscore, 1),
                flag,
            ]
        )

        writer.writerow(row)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert BLAST outfmt 6 alignments into pooled and "
            "length-weighted evolutionary distances."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "blast_table",
        help="BLAST outfmt 6 file, or '-' for stdin",
    )

    parser.add_argument(
        "-o",
        "--out",
        default="-",
        help="output TSV path, or '-' for stdout",
    )

    parser.add_argument(
        "--seqtype",
        choices=["nucl", "prot"],
        default="nucl",
        help="nucleotide (JC69) or protein (Poisson and Kimura) models",
    )

    parser.add_argument(
        "--basis",
        choices=["nogap", "gapped"],
        default="nogap",
        help="p-distance supplied to substitution models",
    )

    parser.add_argument(
        "--aggregate",
        choices=["none", "best", "pool"],
        default="none",
        help=(
            "one row per HSP, best HSP per pair, or pooled per subject locus"
        ),
    )

    parser.add_argument(
        "--locus-gap",
        type=int,
        default=5000,
        help=(
            "maximum subject-coordinate gap merged into one locus when "
            "--aggregate pool is used"
        ),
    )

    parser.add_argument(
        "--min-aln-length",
        type=int,
        default=0,
        help="discard HSPs shorter than this many aligned columns",
    )

    parser.add_argument(
        "--columns",
        default=",".join(DEFAULT_COLUMNS),
        help="comma-separated BLAST outfmt 6 column names in file order",
    )

    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(
        argv if argv is not None else sys.argv[1:]
    )

    if args.locus_gap < 0:
        sys.exit("error: --locus-gap must be >= 0")

    if args.min_aln_length < 0:
        sys.exit("error: --min-aln-length must be >= 0")

    columns = [
        column.strip()
        for column in args.columns.split(",")
        if column.strip()
    ]

    models = (
        NUCL_MODELS
        if args.seqtype == "nucl"
        else PROT_MODELS
    )

    hsps = read_hsps(
        args.blast_table,
        columns,
        args.min_aln_length,
    )

    if not hsps:
        sys.exit("error: no HSPs parsed")

    units = build_units(
        hsps,
        args.aggregate,
        args.locus_gap,
    )

    out_handle = (
        sys.stdout
        if args.out == "-"
        else open(args.out, "w", newline="")
    )

    try:
        write_table(
            units,
            models,
            args.basis,
            out_handle,
        )
    finally:
        if out_handle is not sys.stdout:
            out_handle.close()

    sys.stderr.write(
        f"parsed {len(hsps)} HSPs, wrote {len(units)} rows "
        f"[seqtype={args.seqtype}, basis={args.basis}, "
        f"aggregate={args.aggregate}]\n"
    )


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        try:
            sys.stdout.close()
        except Exception:
            pass
        sys.exit(0)
