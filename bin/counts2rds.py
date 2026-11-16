#!/usr/bin/env python3
"""Merge vcf2counts per-sample CSVs into one analysis-ready RDS.

--input-dir holds one CSV per sample (first column = variant id, the rest one count
column per caller). Each row is tagged with sample_name from its filename, joined
to sample demographics from --metadata, deduplicated on (sample_name, variant_called),
and written as the RDS the diversity analysis expects.

Usage:
    counts2rds.py --input-dir ./vcf2counts/ --metadata sample_metadata.csv \
        --population populations.tsv --output all_samples.rds

--metadata is required and is either analysis-shaped (one row per sample) or the
1000G sequence index (many rows per sample, tab-separated). Population columns
are optional: if absent, metadata supplies sample mapping and custom analysis
columns while population fields default to NA.
Output columns: sample_name, variant_called, <one per caller>, population_code,
super_population_code, population_description, super_population_description,
plus any additional custom columns present in the metadata.

Genotypes: 0 = hom ref, 1 = het, 2 = hom alt, NA = not called (absent from that
sample's CSV, or uncalled/multi-allelic). NA is never filled with 0: that would
claim a homozygous reference call the data never reported.
"""

import argparse
import glob
import gzip
import os
import shutil

import pandas as pd
import pyreadr

# Source column name -> analysis name, first match wins. This is what reads both
# the analysis-shaped metadata and the 1000G sequence index without a flag. Only
# standard population/sample columns that reach the RDS are mapped. Any custom
# analysis columns present in the metadata will be preserved through the merge.
ALIASES = [
    ("sample_name", ["Run", "SAMPLE_NAME", "sample_name", "Sample", "sample"]),
    ("population_code", ["population", "POPULATION", "population_code", "Population"]),
    ("super_population_code", ["Super_Population_Code", "super_population_code",
                               "SUPER_POPULATION_CODE", "super_population"]),
    ("population_description", ["Population_Description", "population_description",
                                "POPULATION_DESCRIPTION"]),
    ("super_population_description", ["Super_Population_Description",
                                      "super_population_description",
                                      "SUPER_POPULATION_DESCRIPTION"]),
]
# The lookup additionally accepts the space-and-title spelling populations.tsv uses.
POP_LOOKUP_ALIASES = [
    (target, names + [n.replace("_", " ").title() for n in names
                      if n.replace("_", " ").title() not in names])
    for target, names in ALIASES
]

ANALYSIS_METADATA_COLS = ["sample_name", "population_code", "super_population_code",
                          "population_description", "super_population_description"]
ID_COLS = ("sample_name", "variant_called")          # not per-caller counts
BLANKS = ("", "na", "nan", "none", "null", ".")      # strings meaning "no value"
SUPER_POP_NAMES = {"AFR": "African", "AMR": "Admixed American", "EAS": "East Asian",
                   "EUR": "European", "SAS": "South Asian"}
_ALIAS_NAMES = {n for _, names in ALIASES + POP_LOOKUP_ALIASES for n in names}


def read_table(path):
    """Read a metadata table.

    The separator comes from the header line, because the sequence index is
    tab-separated with no .tsv extension. If known standard aliases match, we can
    use usecols to optimize, but to fully preserve custom user analysis columns,
    we read the full table structure.
    """
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
    sep = "\t" if "\t" in header else ";" if header.count(";") > header.count(",") else ","
    return pd.read_csv(path, sep=sep, dtype=str, low_memory=False)


def to_analysis(frame, aliases, required, what):
    """Map a frame's columns onto analysis names; abort if one required is absent."""
    mapping = {name: target for target, names in aliases
               for name in names if name in frame.columns}
    out = frame.rename(columns=mapping)
    out = out.loc[:, ~out.columns.duplicated()]   # "population" + "population_code"
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise SystemExit(
            f"{what} is missing required column(s): {', '.join(missing)}\n"
            f"  looking for any of: "
            + "; ".join(f"{t} = {names}" for t, names in aliases if t in required)
            + f"\n  file columns: {', '.join(frame.columns)}")
    return out


def blank_to_na(series):
    """Empty / NA-like strings -> a real NA."""
    return series.astype("object").where(
        ~series.astype(str).str.strip().str.lower().isin(BLANKS), other=pd.NA)


def normalise_metadata(demo, lookup=None):
    """Build analysis metadata columns, preserving any custom user columns.

    The sequence index carries many rows per sample, so it is collapsed to one row
    per sample first: that is what keeps the merge at sample count, not row count.
    """
    demo = to_analysis(demo, ALIASES, ["sample_name"], "Metadata file")

    # Ensure all standard analysis metadata columns exist (population fields default to NA if absent)
    for col in ANALYSIS_METADATA_COLS:
        if col not in demo.columns:
            demo[col] = pd.NA
        demo[col] = blank_to_na(demo[col])

    demo = demo[demo["sample_name"].notna()]
    if demo.empty:
        raise SystemExit("Metadata file has no rows with a sample name")
    n_rows_in = len(demo)

    # Collapse duplicates, keeping the first row per sample
    demo = demo.drop_duplicates(subset="sample_name", keep="first")
    print(f"Metadata: {n_rows_in} rows -> {len(demo)} sample(s) after collapsing")

    has_pop = demo["population_code"].notna().any()

    if lookup is not None and has_pop:
        lookup = to_analysis(lookup, POP_LOOKUP_ALIASES, ["population_code"],
                             "Population lookup file")
        lookup = lookup.drop_duplicates(subset="population_code", keep="first")
        lookup = lookup.set_index("population_code")
        for col in ("super_population_code", "population_description",
                    "super_population_description"):
            if col in lookup.columns:
                demo[col] = demo[col].fillna(demo["population_code"].map(lookup[col]))
        print(f"Population lookup: {len(lookup)} population(s) applied")
    elif lookup is not None and not has_pop:
        print("  NOTE: --population-lookup provided, but metadata has no population codes; lookup ignored")

    if has_pop:
        # populations.tsv carries only the superpopulation CODE, so derive its label.
        demo["super_population_description"] = demo["super_population_description"].fillna(
            demo["super_population_code"].map(SUPER_POP_NAMES))
        if n_filled := int(demo["super_population_description"].notna().sum()):
            print(f"Superpopulation descriptions: filled {n_filled} sample(s) "
                  f"from the superpopulation code")

        # Fill leftovers from another row sharing the population_code.
        for col in ("super_population_code", "super_population_description",
                    "population_description"):
            known = (demo[demo[col].notna()].drop_duplicates(subset="population_code")
                     .set_index("population_code")[col])
            demo[col] = demo[col].fillna(demo["population_code"].map(known))

    return demo.reset_index(drop=True)


def load_sample(path, empty_samples):
    """Read one per-sample CSV, or record it as empty and return None.

    A header-only CSV is a legitimate output (no variant survived) and must not
    abort the run, nor vanish silently: the sample count would then disagree with
    the samplesheet with nothing in the log to explain it.
    """
    sample_name = os.path.splitext(os.path.basename(path))[0]
    df = pd.read_csv(path)
    if df.shape[0] == 0:
        empty_samples.append({"sample_name": sample_name,
                              "source_file": os.path.basename(path),
                              "reason": "no data rows (header only)"})
        return None
    if df.columns[0] != "variant_called":
        df = df.rename(columns={df.columns[0]: "variant_called"})
    df.insert(0, "sample_name", sample_name)
    return df


def save_rds(data, path):
    """Write an RDS whose bytes do not depend on the run.

    pyreadr's gzip writes the mtime and the file name into the gzip header, so the
    same data gives different bytes each run and the nf-test md5 snapshots could
    never match. Recompressing with mtime=0 and an empty filename (what R's saveRDS
    writes) fixes that. Level 6 is the knee of the curve: level 9 is ~4x slower for
    ~3% smaller output.
    """
    output_dir = os.path.dirname(os.path.abspath(path))
    raw_path = os.path.join(output_dir,
                            os.path.splitext(os.path.basename(path))[0] + ".tmp.rds")
    try:
        pyreadr.write_rds(raw_path, data)
        with open(raw_path, "rb") as raw, open(path, "wb") as out, \
                gzip.GzipFile(filename="", fileobj=out, mode="wb",
                              compresslevel=6, mtime=0) as gz:
            shutil.copyfileobj(raw, gz, length=1024 * 1024)
    finally:
        if os.path.exists(raw_path):
            os.remove(raw_path)


def load_counts(input_dir, staged):
    """Read, coerce and dedupe every per-sample CSV in input_dir.

    Returns (data, count_cols, empty_samples, n_dup_keys).
    """
    csv_files = sorted(glob.glob(os.path.join(input_dir, "*.csv")))
    if not csv_files:
        raise SystemExit(f"No CSV files found in {input_dir}")

    # Nextflow stages metadata and lookup into the same work dir as the count CSVs,
    # so the *.csv glob picks them up and pandas aborts on unrecognized structures.
    # Compared by resolved path: they are staged as symlinks.
    staged_paths = {os.path.realpath(f) for f in staged if f}
    excluded = [f for f in csv_files if os.path.realpath(f) in staged_paths]
    if excluded:
        csv_files = [f for f in csv_files if f not in excluded]
        print(f"Excluded {len(excluded)} staged metadata/lookup CSV(s) from the "
              f"count-table glob: {', '.join(sorted(os.path.basename(f) for f in excluded))}")
    if not csv_files:
        raise SystemExit(f"No count tables left in {input_dir}: every .csv "
                         f"there was staged metadata/lookup")
    print(f"Found {len(csv_files)} per-sample CSV(s) in {input_dir}")

    empty_samples, dfs = [], []
    for csv_file in csv_files:
        df = load_sample(csv_file, empty_samples)
        if df is not None:
            dfs.append(df)
    for s in empty_samples:
        print(f"  WARNING: dropping sample {s['sample_name']} - {s['reason']}")
    if not dfs:
        raise SystemExit("No sample CSVs contained any data rows")

    data = pd.concat(dfs, ignore_index=True)
    print(f"Concatenated: {len(data)} rows x {len(data.columns)} cols "
          f"from {len(dfs)} sample(s)")

    # Every non-id column is a caller count. Coerced unconditionally: pandas 3 reads
    # these text CSVs as its `str` dtype, which is not `object`, so a dtype check
    # would skip the conversion and astype would fail on the strings. Int32, not
    # float: pandas needs float to hold NA, but pyreadr writes float as R numeric,
    # while Int32 keeps R integer and still holds NA.
    count_cols = [c for c in data.columns if c not in ID_COLS]
    if not count_cols:
        raise SystemExit(f"No caller count columns found; only id columns "
                         f"{list(ID_COLS)} are present. Check that vcf2counts "
                         f"wrote more than a variant id.")
    data[count_cols] = data[count_cols].apply(pd.to_numeric, errors="coerce").astype("Int32")
    if len(count_cols) == 1:
        print(f"Only one caller column ({count_cols[0]}); nothing to compare "
              f"between callers, writing the single-arm count matrix")

    # Collapse repeated (sample_name, variant_called) by taking the max of each
    # count column: the strongest call wins (2 beats 1 beats 0). Required, because
    # data_io.R builds rownames as paste0(sample_name, "_", variant_called) and
    # data.table rejects duplicates. sort=False keeps the original row order.
    n_dup_keys = int(data.duplicated(subset=list(ID_COLS)).sum())
    if n_dup_keys > 0:
        data = data.groupby(list(ID_COLS), as_index=False, sort=False)[count_cols].max()
        print(f"Deduplicated: collapsed {n_dup_keys} duplicate "
              f"(sample_name, variant_called) key(s); {len(data)} rows remain")
    else:
        print("Deduplicated: no duplicate (sample_name, variant_called) keys found")
    return data, count_cols, empty_samples, n_dup_keys


def attach_metadata(data, metadata, lookup_path):
    """Left-join sample demographics and custom columns onto the counts.

    Returns (data, n_superpop, n_unmatched).
    """
    for path, what in ((metadata, "Metadata file"),
                       (lookup_path, "Population lookup file")):
        if path and not os.path.exists(path):
            raise SystemExit(f"{what} not found: {path}")
    lookup = read_table(lookup_path) if lookup_path else None
    demo = normalise_metadata(read_table(metadata), lookup=lookup)

    has_pop = "population_code" in demo.columns and demo["population_code"].notna().any()
    n_superpop = demo["super_population_code"].nunique(dropna=True) if "super_population_code" in demo.columns else 0

    if has_pop:
        print(f"Metadata: {demo['population_code'].nunique()} population(s), "
              f"{n_superpop} superpopulation(s)")
        if n_superpop == 0:
            print("  NOTE: no superpopulation values available; per-superpopulation "
                  "stages will be empty. Pass --population to supply them.")
    else:
        print("Metadata: loaded sample names and custom columns (no population information present)")

    data = data.merge(demo, on="sample_name", how="left")

    if has_pop:
        n_unmatched = int(data["population_code"].isna().sum())
        if n_unmatched:
            unmatched = sorted(data.loc[data["population_code"].isna(), "sample_name"].unique())
            print(f"  WARNING: {n_unmatched} row(s) from {len(unmatched)} sample(s) "
                  f"have no metadata match: {', '.join(unmatched[:10])}"
                  f"{' ...' if len(unmatched) > 10 else ''}")
            print("           these samples will drop out of the per-population stages")
    else:
        n_matched = int(data["sample_name"].isin(demo["sample_name"]).sum())
        n_unmatched = len(data) - n_matched
        if n_unmatched:
            print(f"  NOTE: {n_unmatched} row(s) have no metadata match.")

    return data, n_superpop, n_unmatched


def save_outputs(data, args, count_cols, empty_samples, n_dup_keys, n_na_cells,
                 n_superpop, n_unmatched):
    """Write the RDS plus the dropped-sample and summary reports."""
    output_dir = os.path.dirname(os.path.abspath(args.output))
    base = os.path.splitext(os.path.basename(args.output))[0]

    # Convert standard metadata columns to nullable string dtype for stable R types
    for col in ANALYSIS_METADATA_COLS:
        if col in data.columns:
            data[col] = data[col].astype("string")

    # pyreadr rejects duplicate column names, and the merge can produce one if the
    # metadata carries a column already present in the counts.
    duplicated = sorted(set(data.columns[data.columns.duplicated()]))
    if duplicated:
        raise SystemExit(f"Duplicate column name(s) in the merged table: "
                         f"{', '.join(duplicated)}. Rename or drop them before "
                         f"writing the RDS.")

    save_rds(data, args.output)
    print(f"Saved RDS: {args.output}")
    if args.write_csv:
        csv_out = os.path.join(output_dir, base + ".csv")
        data.to_csv(csv_out, index=False, na_rep="NA")
        print(f"Saved CSV: {csv_out}")

    # Written even when empty, so the file always exists for MultiQC reporting.
    dropped_path = os.path.join(output_dir, base + "_dropped_samples.csv")
    pd.DataFrame(empty_samples,
                 columns=["sample_name", "source_file", "reason"]).to_csv(
        dropped_path, index=False)
    print(f"Saved dropped-sample report: {dropped_path}"
          f" ({len(empty_samples)} sample(s))")

    n_pops = data["population_code"].nunique() if ("population_code" in data.columns and data["population_code"].notna().any()) else 0

    summary_path = os.path.join(output_dir, base + "_summary.csv")
    pd.DataFrame({
        "metric": ["n_samples", "n_callers", "callers", "n_rows", "n_variants",
                   "n_duplicate_keys_collapsed", "n_samples_dropped",
                   "n_rows_without_metadata", "n_na_genotypes", "n_populations",
                   "n_superpopulations", "metadata_file", "population_file"],
        "value": [data["sample_name"].nunique(), len(count_cols),
                  "|".join(count_cols), len(data), data["variant_called"].nunique(),
                  n_dup_keys, len(empty_samples), n_unmatched, n_na_cells,
                  n_pops, n_superpop, args.metadata,
                  args.population or "none"],
    }).to_csv(summary_path, index=False)
    print(f"Saved summary: {summary_path}")
    print(f"  Rows: {len(data)} | Samples: {data['sample_name'].nunique()} | "
          f"Callers: {count_cols}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-dir", required=True,
                   help="Dir with per-sample vcf2counts CSVs")
    p.add_argument("--metadata", required=True,
                   help="Sample metadata file (required: supplies sample mapping "
                        "and any custom analysis columns; population columns are optional)")
    p.add_argument("--population", default=None,
                   help="Table keyed on population code supplying the superpopulation "
                        "code and the two *_description columns, e.g. populations.tsv")
    p.add_argument("--output", required=True, help="Output RDS file path")
    p.add_argument("--write-csv", action="store_true",
                   help="Also write a CSV copy of the merged table (large; the RDS "
                        "carries the same content)")
    args = p.parse_args()

    data, count_cols, empty_samples, n_dup_keys = load_counts(
        args.input_dir, (args.metadata, args.population))

    # Counted after the dedupe, so the number describes the matrix actually written.
    n_na_cells = int(data[count_cols].isna().sum().sum())
    if n_na_cells:
        n_cells = len(data) * len(count_cols)
        print(f"NA genotypes: {n_na_cells} of {n_cells} cell(s) "
              f"({100 * n_na_cells / n_cells:.1f}%) are NA, i.e. not called: "
              f"the variant is absent from that sample's CSV, or the caller "
              f"left it uncalled/multi-allelic")

    data, n_superpop, n_unmatched = attach_metadata(
        data, args.metadata, args.population)

    save_outputs(data, args, count_cols, empty_samples, n_dup_keys, n_na_cells,
                 n_superpop, n_unmatched)


if __name__ == "__main__":
    main()
