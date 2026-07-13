#!/usr/bin/env python3

import argparse
import csv
import math
from collections import defaultdict


def mean(values):
    return sum(values) / len(values) if values else float("nan")


def median(values):
    if not values:
        return float("nan")

    vals = sorted(values)
    n = len(vals)
    mid = n // 2

    if n % 2:
        return vals[mid]

    return (vals[mid - 1] + vals[mid]) / 2


def sample_sd(values):
    if len(values) < 2:
        return float("nan")

    m = mean(values)

    return math.sqrt(
        sum((x - m) ** 2 for x in values) / (len(values) - 1)
    )


def fmt(value):
    if isinstance(value, float) and math.isnan(value):
        return "NA"

    return f"{value:.6f}"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize BLAST-derived divergence metrics, optionally "
            "by genomic region."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Per-HSP divergence TSV from blast_divergence.py",
    )

    parser.add_argument(
        "-r",
        "--regions",
        default=None,
        help=(
            "Optional region TSV with columns: region, chrom, start, end. "
            "If omitted, all HSPs are summarized together."
        ),
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output summary TSV",
    )

    return parser.parse_args()


def read_regions(path):
    """
    Read region definitions.

    Expected columns:
        region
        chrom
        start
        end
    """

    regions = []

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


def assign_region(chrom, hsp_start, hsp_end, regions):
    """
    Assign an HSP to the first overlapping region.

    If no region file was supplied, all HSPs are assigned to 'all'.
    """

    if not regions:
        return "all"

    for region in regions:
        overlaps = (
            chrom == region["chrom"]
            and hsp_end >= region["start"]
            and hsp_start <= region["end"]
        )

        if overlaps:
            return region["region"]

    return None


def main():
    args = parse_args()

    regions = read_regions(args.regions)

    stats = defaultdict(
        lambda: {
            "d_p_distance": [],
            "d_jc69": [],
            "pident": [],
            "aln_length": [],
            "identities": 0,
            "mismatch": 0,
            "gaps": 0,
        }
    )

    with open(args.input, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_input = {
            "sseqid",
            "s_start",
            "s_end",
            "d_p_distance",
            "d_jc69",
            "pident",
            "aln_length",
            "identities",
            "mismatch",
            "gaps",
        }

        missing_input = required_input - set(reader.fieldnames or [])

        if missing_input:
            raise ValueError(
                "Input divergence file is missing columns: "
                + ", ".join(sorted(missing_input))
            )

        for row in reader:
            chrom = row["sseqid"]

            hsp_start = min(
                int(row["s_start"]),
                int(row["s_end"]),
            )

            hsp_end = max(
                int(row["s_start"]),
                int(row["s_end"]),
            )

            assigned_region = assign_region(
                chrom=chrom,
                hsp_start=hsp_start,
                hsp_end=hsp_end,
                regions=regions,
            )

            # When regions are supplied, skip HSPs outside all regions.
            if assigned_region is None:
                continue

            d_p = row["d_p_distance"]
            d_jc = row["d_jc69"]

            if d_p not in {"", "NA"}:
                stats[assigned_region]["d_p_distance"].append(
                    float(d_p)
                )

            if d_jc not in {"", "NA"}:
                stats[assigned_region]["d_jc69"].append(
                    float(d_jc)
                )

            stats[assigned_region]["pident"].append(
                float(row["pident"])
            )

            stats[assigned_region]["aln_length"].append(
                int(row["aln_length"])
            )

            stats[assigned_region]["identities"] += int(
                row["identities"]
            )

            stats[assigned_region]["mismatch"] += int(
                row["mismatch"]
            )

            stats[assigned_region]["gaps"] += int(
                row["gaps"]
            )

    if not stats:
        raise ValueError(
            "No HSPs were summarized. Check that the input contains data "
            "and, if using --regions, that chromosome names and coordinates "
            "overlap."
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
                "n_hsp",
                "total_aligned_bp",
                "mean_pident",
                "pooled_p_distance_nogap",
                "mean_d_p_distance",
                "median_d_p_distance",
                "sd_d_p_distance",
                "mean_d_jc69",
                "median_d_jc69",
                "sd_d_jc69",
                "total_gaps",
            ]
        )

        for region_name in sorted(stats):
            s = stats[region_name]

            ungapped = (
                s["identities"]
                + s["mismatch"]
            )

            pooled_p = (
                s["mismatch"] / ungapped
                if ungapped > 0
                else float("nan")
            )

            writer.writerow(
                [
                    region_name,
                    len(s["aln_length"]),
                    sum(s["aln_length"]),
                    fmt(mean(s["pident"])),
                    fmt(pooled_p),
                    fmt(mean(s["d_p_distance"])),
                    fmt(median(s["d_p_distance"])),
                    fmt(sample_sd(s["d_p_distance"])),
                    fmt(mean(s["d_jc69"])),
                    fmt(median(s["d_jc69"])),
                    fmt(sample_sd(s["d_jc69"])),
                    s["gaps"],
                ]
            )


if __name__ == "__main__":
    main()
