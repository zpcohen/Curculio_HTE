#!/usr/bin/env python3
"""
curate_divergence_weighted.py

Summarize output from blast_divergence_weighted.py, optionally by genomic
region.

The primary regional estimate is:

    pooled_d_jc69_nogap

This value is calculated by pooling identities and mismatches across every
included alignment row, calculating the site-weighted ungapped p-distance,

    pooled_p = total_mismatches / (total_identities + total_mismatches),

and then applying the Jukes-Cantor 1969 correction once to that pooled
p-distance.

This is preferable to averaging JC69 values across HSPs because every aligned
ungapped nucleotide contributes once to the regional estimate.

The script also reports:
    - the ungapped-length-weighted mean of row-level JC69 distances;
    - unweighted mean, median, and sample SD of row-level distances;
    - pooled percent identity;
    - total aligned and ungapped bases.

Input is expected to be the TSV produced by blast_divergence_weighted.py.
"""

import argparse
import csv
import math
from collections import defaultdict
from typing import Dict, List, Optional


NAN = float("nan")


def mean(values: List[float]) -> float:
    return sum(values) / len(values) if values else NAN


def median(values: List[float]) -> float:
    if not values:
        return NAN

    vals = sorted(values)
    n = len(vals)
    mid = n // 2

    if n % 2:
        return vals[mid]

    return (vals[mid - 1] + vals[mid]) / 2.0


def sample_sd(values: List[float]) -> float:
    if len(values) < 2:
        return NAN

    m = mean(values)
    return math.sqrt(
        sum((value - m) ** 2 for value in values)
        / (len(values) - 1)
    )


def weighted_mean(values: List[float], weights: List[int]) -> float:
    if not values or not weights or len(values) != len(weights):
        return NAN

    denominator = sum(weights)
    if denominator <= 0:
        return NAN

    return sum(
        value * weight
        for value, weight in zip(values, weights)
    ) / denominator


def jc69_distance(p_distance: float) -> float:
    """
    Apply the Jukes-Cantor 1969 correction.

    JC69 is undefined when p >= 0.75 because the logarithm argument is
    zero or negative.
    """
    if math.isnan(p_distance) or p_distance < 0.0:
        return NAN

    argument = 1.0 - (4.0 / 3.0) * p_distance
    if argument <= 0.0:
        return NAN

    return -0.75 * math.log(argument)


def fmt(value: float) -> str:
    if isinstance(value, float) and math.isnan(value):
        return "NA"

    return f"{value:.6f}"


def parse_optional_float(value: str) -> Optional[float]:
    if value in {"", "NA", "nan", "NaN"}:
        return None

    return float(value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize pooled and length-weighted JC69 divergence "
            "from blast_divergence_weighted.py, optionally by region."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="TSV produced by blast_divergence_weighted.py",
    )

    parser.add_argument(
        "-r",
        "--regions",
        default=None,
        help=(
            "Optional region TSV with columns: region, chrom, start, end. "
            "If omitted, all alignment rows are summarized together."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output summary TSV",
    )

    return parser.parse_args()


def read_regions(path: Optional[str]) -> List[Dict[str, object]]:
    """
    Read region definitions.

    Required columns:
        region
        chrom
        start
        end
    """
    regions: List[Dict[str, object]] = []

    if path is None:
        return regions

    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"region", "chrom", "start", "end"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                "Region file is missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:
            start = int(row["start"])
            end = int(row["end"])

            regions.append(
                {
                    "region": row["region"],
                    "chrom": row["chrom"],
                    "start": min(start, end),
                    "end": max(start, end),
                }
            )

    return regions


def assign_region(
    chrom: str,
    aln_start: int,
    aln_end: int,
    regions: List[Dict[str, object]],
) -> Optional[str]:
    """
    Assign an alignment row to the first overlapping region.

    If no region file was supplied, every row is assigned to ``all``.
    """
    if not regions:
        return "all"

    for region in regions:
        overlaps = (
            chrom == region["chrom"]
            and aln_end >= region["start"]
            and aln_start <= region["end"]
        )

        if overlaps:
            return str(region["region"])

    return None


def main() -> None:
    args = parse_args()
    regions = read_regions(args.regions)

    stats = defaultdict(
        lambda: {
            "row_d_p_distance": [],
            "row_d_jc69": [],
            "row_d_jc69_weights": [],
            "row_weighted_mean_d_jc69_hsp": [],
            "pident": [],
            "aln_length": [],
            "ungapped_length": [],
            "identities": 0,
            "mismatch": 0,
            "gaps": 0,
            "n_hsp": 0,
            "n_rows": 0,
        }
    )

    with open(args.input, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_input = {
            "sseqid",
            "s_start",
            "s_end",
            "n_hsp",
            "aln_length",
            "ungapped_length",
            "identities",
            "mismatch",
            "gaps",
            "pident",
            "d_p_distance",
            "d_jc69",
            "weighted_mean_d_jc69_hsp",
        }

        missing_input = required_input - set(reader.fieldnames or [])

        if missing_input:
            raise ValueError(
                "Input divergence file is missing columns: "
                + ", ".join(sorted(missing_input))
            )

        for row in reader:
            chrom = row["sseqid"]

            aln_start = min(
                int(row["s_start"]),
                int(row["s_end"]),
            )
            aln_end = max(
                int(row["s_start"]),
                int(row["s_end"]),
            )

            assigned_region = assign_region(
                chrom=chrom,
                aln_start=aln_start,
                aln_end=aln_end,
                regions=regions,
            )

            # When regions are supplied, skip rows outside all regions.
            if assigned_region is None:
                continue

            s = stats[assigned_region]

            identities = int(row["identities"])
            mismatches = int(row["mismatch"])
            gaps = int(row["gaps"])
            aln_length = int(row["aln_length"])
            ungapped_length = int(row["ungapped_length"])
            n_hsp = int(row["n_hsp"])

            d_p = parse_optional_float(row["d_p_distance"])
            d_jc69 = parse_optional_float(row["d_jc69"])
            weighted_hsp_jc69 = parse_optional_float(
                row["weighted_mean_d_jc69_hsp"]
            )

            if d_p is not None:
                s["row_d_p_distance"].append(d_p)

            if d_jc69 is not None and ungapped_length > 0:
                s["row_d_jc69"].append(d_jc69)
                s["row_d_jc69_weights"].append(ungapped_length)

            if weighted_hsp_jc69 is not None:
                s["row_weighted_mean_d_jc69_hsp"].append(
                    weighted_hsp_jc69
                )

            s["pident"].append(float(row["pident"]))
            s["aln_length"].append(aln_length)
            s["ungapped_length"].append(ungapped_length)
            s["identities"] += identities
            s["mismatch"] += mismatches
            s["gaps"] += gaps
            s["n_hsp"] += n_hsp
            s["n_rows"] += 1

    if not stats:
        raise ValueError(
            "No alignment rows were summarized. Check that the input "
            "contains data and, if using --regions, that chromosome names "
            "and coordinates overlap."
        )

    with open(args.output, "w", newline="") as out:
        writer = csv.writer(
            out,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "region",
                "n_rows",
                "n_hsp",
                "total_aligned_bp",
                "total_ungapped_bp",
                "pooled_pident_nogap",
                "pooled_p_distance_nogap",
                "pooled_d_jc69_nogap",
                "weighted_mean_d_jc69_by_row",
                "mean_d_p_distance_by_row",
                "median_d_p_distance_by_row",
                "sd_d_p_distance_by_row",
                "mean_d_jc69_by_row",
                "median_d_jc69_by_row",
                "sd_d_jc69_by_row",
                "mean_weighted_mean_d_jc69_hsp_by_row",
                "total_identities",
                "total_mismatches",
                "total_gaps",
                "jc69_flag",
            ]
        )

        for region_name in sorted(stats):
            s = stats[region_name]

            ungapped = s["identities"] + s["mismatch"]

            pooled_p = (
                s["mismatch"] / ungapped
                if ungapped > 0
                else NAN
            )

            pooled_pident = (
                100.0 * s["identities"] / ungapped
                if ungapped > 0
                else NAN
            )

            pooled_jc69 = jc69_distance(pooled_p)

            row_weighted_jc69 = weighted_mean(
                s["row_d_jc69"],
                s["row_d_jc69_weights"],
            )

            jc69_flag = (
                "saturated"
                if isinstance(pooled_jc69, float)
                and math.isnan(pooled_jc69)
                and not math.isnan(pooled_p)
                else ""
            )

            writer.writerow(
                [
                    region_name,
                    s["n_rows"],
                    s["n_hsp"],
                    sum(s["aln_length"]),
                    sum(s["ungapped_length"]),
                    fmt(pooled_pident),
                    fmt(pooled_p),
                    fmt(pooled_jc69),
                    fmt(row_weighted_jc69),
                    fmt(mean(s["row_d_p_distance"])),
                    fmt(median(s["row_d_p_distance"])),
                    fmt(sample_sd(s["row_d_p_distance"])),
                    fmt(mean(s["row_d_jc69"])),
                    fmt(median(s["row_d_jc69"])),
                    fmt(sample_sd(s["row_d_jc69"])),
                    fmt(mean(s["row_weighted_mean_d_jc69_hsp"])),
                    s["identities"],
                    s["mismatch"],
                    s["gaps"],
                    jc69_flag,
                ]
            )


if __name__ == "__main__":
    main()
