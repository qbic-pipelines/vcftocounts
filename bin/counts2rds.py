#!/usr/bin/env python3
"""Merge vcf2counts per-sample CSVs into one analysis-ready RDS.

Reads all *.csv from --input-dir (one per sample, first column = variant id,
remaining columns = one count column per caller), adds sample_name from the
filename, concatenates, joins sample demographics, deduplicates
(sample_name, variant_called) pairs, and saves as RDS.

The RDS is intended to be ready for the diversity analysis: the demographic
columns carry the population_code / super_population_code / *_description
naming that the R functions expect, so no further renaming is needed.

Usage:
    counts2rds.py --input-dir ./vcf2counts/ --metadata demo.csv --output all_samples.rds

Input (per-sample CSV):
    variant,sarek_freebayes,sarek_strelka
    chr22_00001,1,0
    chr22_00002,1,1

--metadata accepts either of two shapes, detected from the column names:

  1. Analysis-shaped file (one row per sample):
       Run,population,Super_Population_Code,Population_Description,Super_Population_Description

  2. 1000 Genomes sequence index (many rows per sample, tab-separated):
       FASTQ_FILE  MD5  RUN_ID  ...  SAMPLE_NAME  POPULATION  ...
     Only SAMPLE_NAME and POPULATION are read; the other 24 columns are ignored,
     so passing the 67 MB index costs one two-column scan.

Superpopulation values are not present in the sequence index. Supply them with
--population-lookup, a small table keyed on population code:

    population,super_population_code,population_description,super_population_description
    GBR,EUR,British,European

Output (one RDS):
    sample_name,variant_called,sarek_freebayes,sarek_strelka,
    population_code,super_population_code,population_description,super_population_description
    HG00097,chr22_00001,1,0,GBR,EUR,British,European
"""

import argparse
import glob
import os

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import r


# Accepted source column names -> analysis-side names. First match wins, so the
# more specific names are listed before the generic ones. This is what lets one
# script read both the analysis-shaped metadata file and the 1000G sequence
# index without a flag.
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

# Population-code lookup: same idea, for --population-lookup files.
POP_LOOKUP_ALIASES = [
    ("population_code", ["population", "POPULATION", "population_code", "Population"]),
    ("super_population_code", ["super_population_code", "Super_Population_Code",
                               "SUPER_POPULATION_CODE", "super_population"]),
    ("population_description", ["population_description", "Population_Description",
                                "POPULATION_DESCRIPTION"]),
    ("super_population_description", ["super_population_description",
                                      "Super_Population_Description",
                                      "SUPER_POPULATION_DESCRIPTION"]),
]

# The five columns the downstream R code understands, in output order.
ANALYSIS_METADATA_COLS = ["sample_name", "population_code", "super_population_code",
                          "population_description", "super_population_description"]

# Columns of the merged count table that are not per-caller counts.
ID_COLS = ("sample_name", "variant_called")

# Strings that mean "no value" in a metadata file.
BLANKS = ("", "na", "nan", "none", "null", ".")

# Only these columns are ever read from a metadata file. Everything else in a
# sequence index (STUDY_NAME, ANALYSIS_GROUP, CENTER_NAME, ...) is skipped at
# parse time rather than loaded and dropped, which is what keeps the 67 MB
# index cheap.
_ALL_ALIAS_NAMES = {n for _, names in ALIASES for n in names} | \
                   {n for _, names in POP_LOOKUP_ALIASES for n in names}


def detect_separator(path):
    """Return the field separator of a delimited text file.

    The 1000G sequence index is tab-separated but has no .tsv extension, so the
    separator is decided from the header line rather than the filename.
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
    if "\t" in header:
        return "\t"
    if header.count(",") >= header.count(";") and "," in header:
        return ","
    if ";" in header:
        return ";"
    return "\t" if "\t" in header else ","


def read_delimited(path, nrows=None):
    """Read a delimited metadata file, loading only columns we could use.

    Two-pass: read the header, keep the intersection with the alias sets, then
    parse only those columns. For the sequence index this cuts the read from 26
    columns to 2.
    """
    sep = detect_separator(path)
    header = pd.read_csv(path, sep=sep, nrows=0, dtype=str)
    wanted = [c for c in header.columns if c in _ALL_ALIAS_NAMES]
    # If nothing matched, fall back to a full read so the error message can list
    # what the file actually contains.
    usecols = wanted if wanted else None
    return pd.read_csv(path, sep=sep, usecols=usecols, nrows=nrows,
                       dtype=str, low_memory=False)


def rename_to_analysis(frame, aliases, required, what):
    """Map a frame's columns onto analysis names; error if a required one is absent."""
    mapping = {}
    for target, names in aliases:
        for name in names:
            if name in frame.columns:
                mapping[name] = target
                break

    out = frame.rename(columns=mapping)
    # A file listing both "population" and "population_code" would collide; keep
    # the first occurrence.
    out = out.loc[:, ~out.columns.duplicated()]

    missing = [c for c in required if c not in out.columns]
    if missing:
        raise SystemExit(
            f"{what} is missing required column(s): {', '.join(missing)}\n"
            f"  looking for any of: "
            + "; ".join(f"{t} = {names}" for t, names in aliases if t in required)
            + f"\n  file columns: {', '.join(frame.columns)}"
        )
    return out


def blank_to_na(series):
    """Map empty / NA-like strings to a real NA."""
    cleaned = series.astype("object").where(
        ~series.astype(str).str.strip().str.lower().isin(BLANKS), other=pd.NA
    )
    return cleaned


def normalise_metadata(demo, lookup=None):
    """Build the five analysis metadata columns, one row per sample.

    Handles both supported metadata shapes. The sequence index carries many rows
    per sample, so it is collapsed to one row per sample before any join - that
    is what keeps the downstream merge at sample count instead of row count.
    """
    demo = rename_to_analysis(demo, ALIASES,
                              required=["sample_name", "population_code"],
                              what="Metadata file")

    for col in ANALYSIS_METADATA_COLS:
        demo[col] = blank_to_na(demo[col]) if col in demo.columns else pd.NA

    # Drop rows that name no sample or no population before collapsing.
    demo = demo[demo["sample_name"].notna() & demo["population_code"].notna()]
    if demo.empty:
        raise SystemExit("Metadata file has no rows with both a sample and a population")

    n_rows_in = len(demo)

    # Collapse to one row per sample. sort=False keeps first-seen order, which
    # for a sequence index means the first FASTQ entry for that sample.
    demo = demo.drop_duplicates(subset="sample_name", keep="first")
    print(f"Metadata: {n_rows_in} rows -> {len(demo)} sample(s) after collapsing")

    # Fill superpopulation and descriptions from a population-code lookup, if given.
    if lookup is not None:
        lookup = rename_to_analysis(lookup, POP_LOOKUP_ALIASES,
                                    required=["population_code"],
                                    what="Population lookup file")
        lookup = lookup.drop_duplicates(subset="population_code", keep="first")
        lookup = lookup.set_index("population_code")

        for col in ("super_population_code", "population_description",
                    "super_population_description"):
            if col in lookup.columns:
                # .map on the deduplicated column: ~one lookup per population,
                # not per sample row.
                mapped = demo["population_code"].map(lookup[col])
                demo[col] = demo[col].fillna(mapped)
        print(f"Population lookup: {len(lookup)} population(s) applied")

    # Last resort for the analysis-shaped file: borrow a missing superpopulation
    # or description from another row sharing the same population_code.
    incomplete = demo["super_population_code"].isna() | demo["super_population_description"].isna()
    if incomplete.any():
        lookup = (demo[~incomplete]
                  .drop_duplicates(subset="population_code", keep="first")
                  .set_index("population_code")[["super_population_code",
                                                "super_population_description"]])
        for col in ("super_population_code", "super_population_description"):
            demo[col] = demo[col].fillna(demo["population_code"].map(lookup[col]))

    missing_desc = demo["population_description"].isna()
    if missing_desc.any():
        lookup = (demo[~missing_desc]
                  .drop_duplicates(subset="population_code", keep="first")
                  .set_index("population_code")["population_description"])
        demo["population_description"] = demo["population_description"].fillna(
            demo["population_code"].map(lookup)
        )

    # Exactly the five analysis columns, in order. Anything else in the metadata
    # file is dropped here: hill_calculations.R builds its caller-column list as
    # setdiff(colnames, get_metadata_cols()), so a leftover character column
    # would be treated as a caller and break as.matrix().
    out = demo[ANALYSIS_METADATA_COLS].reset_index(drop=True)

    return out


def load_sample(path, empty_samples):
    """Read one per-sample CSV, or record it as empty and return None.

    A header-only CSV is a legitimate Nextflow output (the sample simply had no
    variants surviving), so it must not abort the run. It must not be silently
    dropped either, because the sample count then differs from the samplesheet
    with nothing in the log to explain it. Recorded and reported instead.
    """
    sample_name = os.path.splitext(os.path.basename(path))[0]
    df = pd.read_csv(path)

    if df.shape[0] == 0:
        empty_samples.append(
            {"sample_name": sample_name, "source_file": os.path.basename(path),
             "reason": "no data rows (header only)"}
        )
        return None

    # First column is the variant id; give it the analysis-side name.
    if df.columns[0] != "variant_called":
        df = df.rename(columns={df.columns[0]: "variant_called"})

    # Prepend sample_name so column order is stable regardless of input.
    df.insert(0, "sample_name", sample_name)
    return df


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input-dir", required=True, help="Dir with per-sample vcf2counts CSVs")
    p.add_argument("--metadata", required=False, default=None,
                   help="Sample metadata: analysis-shaped CSV, or the 1000G "
                        "sequence index (SAMPLE_NAME + POPULATION are read)")
    p.add_argument("--population-lookup", required=False, default=None,
                   help="Optional table keyed on population code supplying "
                        "super_population_code and the two *_description columns")
    p.add_argument("--output", required=True, help="Output RDS file path")
    p.add_argument("--write-csv", action="store_true",
                   help="Also write a CSV copy of the merged table. Off by default: "
                        "the CSV is large and the RDS carries the same content.")
    args = p.parse_args()

    # ---------------------------------------------------------------- load
    csv_files = sorted(glob.glob(os.path.join(args.input_dir, "*.csv")))
    if not csv_files:
        raise SystemExit(f"No CSV files found in {args.input_dir}")

    # The metadata and population-lookup files are staged into the same work
    # directory as the per-sample CSVs (--input-dir '.'), so a plain "*.csv"
    # glob picks them up as if they were count tables. That fails confusingly:
    # the metadata is tab-separated with commas inside its description fields,
    # so pandas aborts with "Expected 11 fields in line 258, saw 12" on a file
    # that was never a count table. Excluded by resolved path, because Nextflow
    # stages them as symlinks and the literal argument would not match.
    staged_inputs = {os.path.realpath(p)
                     for p in (args.metadata, args.population_lookup) if p}
    excluded = [f for f in csv_files if os.path.realpath(f) in staged_inputs]
    if excluded:
        csv_files = [f for f in csv_files if os.path.realpath(f) not in staged_inputs]
        print(f"Excluded {len(excluded)} staged metadata/lookup CSV(s) from the "
              f"count-table glob: {', '.join(sorted(os.path.basename(f) for f in excluded))}")
        if not csv_files:
            raise SystemExit(
                f"Every .csv in {args.input_dir} was excluded as metadata/lookup; "
                f"no count tables remain to read")

    print(f"Found {len(csv_files)} per-sample CSV(s) in {args.input_dir}")

    empty_samples = []
    dfs = []
    for csv_file in csv_files:
        df = load_sample(csv_file, empty_samples)
        if df is not None:
            dfs.append(df)

    if not dfs:
        raise SystemExit("No sample CSVs contained any data rows")
    for s in empty_samples:
        print(f"  WARNING: dropping sample {s['sample_name']} - {s['reason']}")

    data = pd.concat(dfs, ignore_index=True, sort=False)
    print(f"Concatenated: {len(data)} rows x {len(data.columns)} cols "
          f"from {len(dfs)} sample(s)")

    # ----------------------------------------------------- fill + coerce
    # Non-numeric columns other than the ids are coerced to int, so a sample
    # whose CSV came back as strings does not break the merge downstream.
    # NOTE: NA genotypes (multi-allelic, uncalled) are written as 0 here,
    # i.e. they are recorded as homozygous reference.
    for col in data.columns:
        if col not in ID_COLS:
            if data[col].dtype == "object":
                data[col] = pd.to_numeric(data[col], errors="coerce")
            data[col] = data[col].fillna(0).astype(int)

    count_cols = [c for c in data.columns if c not in ID_COLS]
    if len(count_cols) < 2:
        raise SystemExit(
            f"Need at least two caller columns to compare; found {count_cols}"
        )

    # ------------------------------------------------------------ dedupe
    # Collapse repeated (sample_name, variant_called) pairs by taking the max of
    # each count column, i.e. the strongest genotype call wins (2 beats 1 beats 0).
    # Required: data.table rejects duplicate row names, and data_io.R assigns
    # rownames as paste0(sample_name, "_", variant_called).
    n_before = len(data)
    # duplicated() is a single hash pass; calling groupby().ngroups here instead
    # would run the whole grouping twice on the common no-duplicates path.
    n_dup_keys = int(data.duplicated(subset=list(ID_COLS)).sum())
    if n_dup_keys > 0:
        data = (data.groupby(list(ID_COLS), as_index=False, sort=False)[count_cols]
                    .max())
        print(f"Deduplicated: collapsed {n_dup_keys} duplicate "
              f"(sample_name, variant_called) key(s); {len(data)} rows remain")
    else:
        print("Deduplicated: no duplicate (sample_name, variant_called) keys found")

    # ---------------------------------------------------------- metadata
    n_unmatched = 0
    n_superpop = 0
    if args.metadata:
        if not os.path.exists(args.metadata):
            raise SystemExit(f"Metadata file not found: {args.metadata}")

        lookup = None
        if args.population_lookup:
            if not os.path.exists(args.population_lookup):
                raise SystemExit(f"Population lookup file not found: {args.population_lookup}")
            lookup = read_delimited(args.population_lookup)

        demo = normalise_metadata(read_delimited(args.metadata), lookup=lookup)
        n_superpop = demo["super_population_code"].nunique(dropna=True)
        print(f"Metadata: {demo['population_code'].nunique()} population(s), "
              f"{n_superpop} superpopulation(s)")
        if n_superpop == 0:
            print("  NOTE: no superpopulation values available; per-superpopulation "
                  "stages will be empty. Pass --population-lookup to supply them.")

        data = data.merge(demo, on="sample_name", how="left")

        n_unmatched = int(data["population_code"].isna().sum())
        if n_unmatched:
            unmatched = sorted(data.loc[data["population_code"].isna(),
                                        "sample_name"].unique())
            print(f"  WARNING: {n_unmatched} row(s) from {len(unmatched)} sample(s) "
                  f"have no metadata match: {', '.join(unmatched[:10])}"
                  f"{' ...' if len(unmatched) > 10 else ''}")
            print("           these samples will drop out of the per-population stages")
    else:
        print("Metadata: none supplied; the RDS carries count columns only")


    # ------------------------------------------------------------- write
    output_dir = os.path.dirname(os.path.abspath(args.output))
    base = os.path.splitext(os.path.basename(args.output))[0]

    # Write the RDS through a temporary CSV instead of converting the pandas
    # frame with rpy2. Direct conversion costs ~100 s for 3.5 M rows: rpy2 walks
    # every column through pandas2ri and copies the mixed character/numeric frame
    # repeatedly, which also produces "Error while trying to convert the column"
    # warnings for all-NA columns. Writing the CSV and letting R parse it with
    # explicit colClasses takes ~16 s and removes the warnings, because the
    # column types are declared rather than guessed. The CSV is transport only:
    # it is deleted afterwards and is not an output of the module.
    transport_csv = os.path.join(output_dir, base + ".transport.tmp.csv")
    col_classes = ["integer" if str(dt).startswith(("int", "uint")) else "character"
                   for dt in data.dtypes]
    try:
        # na_rep="NA" matches read.csv's default na.strings, so missing metadata
        # arrives in R as a true NA rather than an empty string.
        data.to_csv(transport_csv, index=False, na_rep="NA")

        r.assign("transport_csv", transport_csv)
        r.assign("transport_classes", ro.StrVector(col_classes))
        r.assign("transport_out", args.output)
        # comment.char="" is required because "#" would otherwise truncate any
        # line containing one, silently dropping data.
        r("""
          transport_df <- read.csv(transport_csv, header = TRUE,
                                   stringsAsFactors = FALSE,
                                   colClasses = transport_classes,
                                   comment.char = "", check.names = FALSE)
          saveRDS(transport_df, file = transport_out)
        """)
    finally:
        if os.path.exists(transport_csv):
            os.remove(transport_csv)

    print(f"Saved RDS: {args.output}")
    if args.write_csv:
        csv_out = os.path.join(output_dir, base + ".csv")
        data.to_csv(csv_out, index=False)
        print(f"Saved CSV: {csv_out}")

    # Samples dropped for having no rows. Written even when empty, so the file
    # always exists for downstream / MultiQC reporting.
    dropped_path = os.path.join(output_dir, base + "_dropped_samples.csv")
    pd.DataFrame(empty_samples,
                 columns=["sample_name", "source_file", "reason"]).to_csv(
        dropped_path, index=False)
    print(f"Saved dropped-sample report: {dropped_path}"
          f" ({len(empty_samples)} sample(s))")

    summary_path = os.path.join(output_dir, base + "_summary.csv")
    pd.DataFrame(
        {
            "metric": ["n_samples", "n_callers", "callers", "n_rows",
                       "n_variants", "n_duplicate_keys_collapsed",
                       "n_samples_dropped", "n_rows_without_metadata",
                       "n_populations", "n_superpopulations",
                       "metadata_file", "population_lookup_file"],
            "value": [data["sample_name"].nunique(), len(count_cols),
                      "|".join(count_cols), len(data),
                      data["variant_called"].nunique(), n_dup_keys,
                      len(empty_samples), n_unmatched,
                      data["population_code"].nunique() if args.metadata else 0,
                      n_superpop, args.metadata if args.metadata else "none",
                      args.population_lookup if args.population_lookup else "none"],
        }
    ).to_csv(summary_path, index=False)
    print(f"Saved summary: {summary_path}")

    print(f"  Rows: {len(data)} | Samples: {data['sample_name'].nunique()} | "
          f"Callers: {count_cols}")


if __name__ == "__main__":
    main()
