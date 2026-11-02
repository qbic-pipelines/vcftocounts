#!/usr/bin/env python3
"""
VCF to count matrix converter — streaming, memory-efficient.

Usage:
    vcf2counts.py --help
    vcf2counts.py --output=<ofile> VCF

Options:
    -h, --help                  help screen
    -o, --output=<ofile>        output file name [default: mat.csv]

Arguments:
    VCF                     input VCF file (plain or .gz)

The output is a CSV with the same semantics as the original R vcf2counts:
  - Rows are variants labelled "chrom:pos_REF/ALT" (first ALT only)
  - Columns are samples, one per VCF sample column
  - Values: 0 = hom-ref (00), 1 = het (01/10), 2 = hom-alt (11),
    NA = anything else (uncalled, multi-allelic, etc.)
  - The GT field is stripped of all punctuation (/|) before classification,
    so "0/1" → "01", "1|1" → "11", "0/2" → "02" → NA, etc.
"""

import sys
import gzip
import csv
import re
import argparse

# ---- helpers ----

def open_vcf(path):
    """Open a VCF file, transparently handling .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def classify_gt(gt_raw: str) -> int | None:
    """
    Classify a raw GT string the same way as the original R code.

    The R code does:
        xstrip <- gsub("[[:punct:]]", "", x)

    [[:punct:]] strips /, |, . etc. So:
        "0/0"  → "00"  → 0  (hom-ref)
        "0/1"  → "01"  → 1  (het)
        "1|0"  → "10"  → 1  (het, phased)
        "1/1"  → "11"  → 2  (hom-alt)
        "./."  → ".."  → NA (uncalled)
        "0/2"  → "02"  → NA (multi-allelic)

    Then:
        "11"      → 2
        "01"/"10" → 1
        "00"/""   → 0
        anything else → None (→ NA)
    """
    gt_stripped = re.sub(r'[^\w]', '', gt_raw)  # keep only alphanumeric
    if gt_stripped == "11":
        return 2
    if gt_stripped in ("01", "10"):
        return 1
    if gt_stripped in ("00", ""):
        return 0
    return None  # NA


def process_vcf(vcf_path: str, output_path: str):
    """Stream through the VCF line-by-line and write a CSV row as we go."""

    with open_vcf(vcf_path) as fh, \
         open(output_path, "w", newline="") as out:

        writer = None
        sample_names: list[str] | None = None

        for line in fh:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                # Parse header — samples start at column 10 (index 9)
                fields = line.rstrip("\n").split("\t")
                sample_names = fields[9:]
                writer = csv.writer(out)
                # Write header row: empty first cell, then sample names
                writer.writerow([""] + sample_names)
                continue

            # Data line
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue  # skip malformed

            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]

            # Build row label: "chrom:pos_REF/ALT" using only the first ALT
            first_alt = alt.split(",")[0]
            label = f"{chrom}:{pos}_{ref}/{first_alt}"

            # Parse each sample's GT (first sub-field of each sample column)
            gt_fields = fields[9:]  # FORMAT + sample columns
            row = [label]
            for i, sample in enumerate(sample_names):
                if i < len(gt_fields) and gt_fields[i] != ".":
                    gt_raw = gt_fields[i].split(":")[0]  # take GT part only
                    val = classify_gt(gt_raw)
                    row.append("NA" if val is None else str(val))
                else:
                    row.append("")

            writer.writerow(row)

        # Empty VCF (no #CHROM header): write header-only so downstream always
        # sees a file with a recognizable header row instead of a zero-byte
        # file that pandas would interpret as "no columns".
        if writer is None:
            writer = csv.writer(out)
            writer.writerow([""] + list(sample_names or []))


def main():
    parser = argparse.ArgumentParser(
        description="Convert VCF to genotype count matrix (memory-efficient streaming).",
    )
    parser.add_argument("--output", "-o", default="mat.csv",
                        help="Output CSV file (default: mat.csv)")
    parser.add_argument("vcf", help="Input VCF file (plain text or gzip-compressed)")
    args = parser.parse_args()

    process_vcf(args.vcf, args.output)


if __name__ == "__main__":
    main()