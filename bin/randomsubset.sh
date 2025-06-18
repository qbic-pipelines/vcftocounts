#!/usr/bin/env bash

set -euo pipefail

# Usage: randomsubset.sh <input.vcf> <output.vcf> <fraction>
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <input.vcf> <output.vcf> <fraction (e.g. 0.1 for 10%)>"
    exit 1
fi

input_vcf="$1"
output_vcf="$2"
fraction="$3"

# Create temp files and directories
tmpdir=$(mktemp -d)
tmp_vcf="$tmpdir/tmp.vcf"
tmp_sorted_vcf="$tmpdir/tmp.sorted.vcf"

# Calculate number of records to sample
subset_count=$(bcftools stats "$input_vcf" | awk -v frac="$fraction" -F'\t' '$3=="number of records:" {print int($4*frac)}')

echo "Sampling $subset_count records from $input_vcf"

# Write header
bcftools view --header-only "$input_vcf" > "$tmp_vcf"

# Randomly sample records
bcftools view --no-header "$input_vcf" | \
    awk '{printf("%f\t%s\n",rand(),$0);}' | \
    sort -t $'\t' -T "$tmpdir" -k1,1g | \
    head -n "$subset_count" | \
    cut -f 2- >> "$tmp_vcf" || true

# Sort and write to output
bcftools sort -T "$tmpdir" -o "$output_vcf" "$tmp_vcf"

# Clean up
rm -rf "$tmpdir"
