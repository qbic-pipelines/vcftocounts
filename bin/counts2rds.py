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
    counts2rds.py --input-dir ./vcf2counts/ --metadata sequence.index \
        --population populations.tsv --output all_samples.rds

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
--population (also spelled --population-lookup), a small table keyed on
population code, which is exactly the 1000G populations.tsv:

    Population Description  Population Code  Super Population  ...
    Chinese Dai in Xishuangbanna, China  CDX  EAS  ...

Note the space-separated headers; both space and underscore spellings are
accepted. This file has no superpopulation description column, so
super_population_description is derived from the code (EAS -> East Asian).

Output (one RDS):
    sample_name,variant_called,sarek_freebayes,sarek_strelka,
    population_code,super_population_code,population_description,super_population_description
    HG00097,chr22_00001,1,0,GBR,EUR,British,European
"""

import argparse
import glob
import gzip
import os
import shutil

import pandas as pd
import pyreadr

# Accepted source column names -> analysis-side names. First match wins, so the
# more specific names are listed before the generic ones. This is what lets one
# script read both the analysis-shaped metadata file and the 1000G sequence
# index without a flag.
#
# Only columns that end up in the RDS are mapped. Other identifiers a metadata
# file may carry (id, sample_accession, ...) are left alone: nothing downstream
# reads them, and the R count-matrix code derives the caller columns as
# setdiff(colnames, metadata_cols), so an extra character column would be taken
# for a caller. Overwriting them with a flag can come with the first downstream
# consumer that needs them.
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

# Population-code lookup: same idea, for --population files.
#
# The 1000G populations.tsv spells its headers with spaces AND capitals
# ("Population Code", "Super Population", "Population Description"), so those
# exact forms are listed first; the underscore forms cover hand-written tables.
POP_LOOKUP_ALIASES = [
    ("population_code", ["population", "POPULATION", "population_code", "Population",
                         "Population Code"]),
    ("super_population_code", ["super_population_code", "Super_Population_Code",
                               "SUPER_POPULATION_CODE", "super_population",
                               "Super Population", "Super Population Code"]),
    ("population_description", ["population_description", "Population_Description",
                                "POPULATION_DESCRIPTION", "Population Description"]),
    ("super_population_description", ["super_population_description",
                                      "Super_Population_Description",
                                      "SUPER_POPULATION_DESCRIPTION",
                                      "Super Population Description"]),
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
    if ";" in header and header.count(";") > header.count(","):
        return ";"
    return ","


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

    # Collapse to one row per sample. drop_duplicates keeps first-seen order,
    # which for a sequence index means the first FASTQ entry for that sample.
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

    # The 1000G populations.tsv gives only the superpopulation CODE (EAS, EUR,
    # ...), so super_population_description is still empty at this point. Fill it
    # from any other row carrying the same code, which is how a table that has no
    # description column at all still ends up with one. Tone is left minimal: it
    # is a grouping label for per-superpopulation output, nothing downstream
    # parses it.
    SUPER_POP_NAMES = {
        "AFR": "African", "AMR": "Admixed American", "EAS": "East Asian",
        "EUR": "European", "SAS": "South Asian",
    }
    needs_super_desc = (demo["super_population_code"].notna()
                        & demo["super_population_description"].isna())
    if needs_super_desc.any():
        demo["super_population_description"] = demo[
            "super_population_description"
        ].fillna(demo["super_population_code"].map(SUPER_POP_NAMES))
        n_filled = int(demo["super_population_description"].notna().sum())
        if n_filled:
            print(f"Superpopulation descriptions: filled {n_filled} "
                  f"sample(s) from the superpopulation code")

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
    # --population is the name the Nextflow module uses; --population-lookup is
    # the original. Both are accepted (the module's flag wins if both are given)
    # rather than renaming either, so a pipeline run and a hand-typed command
    # both work.
    p.add_argument("--population-lookup", "--population", dest="population_lookup",
                   required=False, default=None,
                   help="Table keyed on population code supplying the "
                        "superpopulation code and the two *_description columns, "
                        "e.g. the 1000G populations.tsv")
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

    data = pd.concat(dfs, ignore_index=True)
    print(f"Concatenated: {len(data)} rows x {len(data.columns)} cols "
          f"from {len(dfs)} sample(s)")

    # ----------------------------------------------------- fill + coerce
    # Every non-id column is a caller count column and is coerced to int, so a
    # sample whose CSV came back as strings does not break the merge downstream.
    # The conversion runs unconditionally: pandas 3 reads the text CSVs as its
    # `str` dtype, which is not `object`, so a dtype check here would skip the
    # conversion and `astype` would fail on the remaining strings.
    # NOTE: NA genotypes (multi-allelic, uncalled) are written as 0 here,
    # i.e. they are recorded as homozygous reference.
    count_cols = [c for c in data.columns if c not in ID_COLS]
    n_na_genotypes = 0
    for col in count_cols:
        numeric = pd.to_numeric(data[col], errors="coerce")
        n_na_genotypes += int(numeric.isna().sum())
        # int32 rather than the default int64: counts are only 0/1/2, and R has
        # no 64-bit integer, so pyreadr writes int64 as a double. Declaring
        # int32 keeps the RDS column an R integer, which is what the diversity
        # analysis expects (and halves the memory of the merge).
        data[col] = numeric.fillna(0).astype("int32")
    if n_na_genotypes:
        print(f"  WARNING: {n_na_genotypes} NA genotype(s) were written as 0 "
              f"(homozygous reference)")

    if not count_cols:
        raise SystemExit(
            f"No caller count columns found; only id columns {list(ID_COLS)} "
            f"are present. Check that vcf2counts wrote more than a variant id."
        )
    if len(count_cols) == 1:
        # One label per sample is a supported input: the samplesheet allows it
        # and BCFTOOLS_MERGE is skipped for a single VCF per sample. The table is
        # still a valid count matrix, just with nothing to compare against.
        print(f"Only one caller column ({count_cols[0]}); nothing to compare "
              f"between callers, writing the single-arm count matrix")

    # ------------------------------------------------------------ dedupe
    # Collapse repeated (sample_name, variant_called) pairs by taking the max of
    # each count column, i.e. the strongest genotype call wins (2 beats 1 beats 0).
    # Required: data.table rejects duplicate row names, and data_io.R assigns
    # rownames as paste0(sample_name, "_", variant_called).
    # duplicated() is a single hash pass; calling groupby().ngroups here instead
    # would run the whole grouping twice on the common no-duplicates path.
    n_dup_keys = int(data.duplicated(subset=list(ID_COLS)).sum())
    if n_dup_keys > 0:
        # groupby sorts its keys unless told otherwise; sort=False keeps the
        # original row order, so the RDS row order does not depend on how many
        # duplicates happened to be present.
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

    # Metadata columns are declared as pandas' nullable string dtype so that an
    # all-empty column (no superpopulation anywhere) is written as R character
    # NA. An object column that is entirely NA is written as logical instead,
    # which would change the column type of the RDS between runs.
    for col in ANALYSIS_METADATA_COLS:
        if col in data.columns:
            data[col] = data[col].astype("string")

    # pyreadr rejects duplicate column names outright, and the merge above can
    # produce them if a metadata file carries a column already in the counts.
    duplicated = sorted(set(data.columns[data.columns.duplicated()]))
    if duplicated:
        raise SystemExit(
            f"Duplicate column name(s) in the merged table: {', '.join(duplicated)}. "
            f"Rename or drop them before writing the RDS."
        )

    # The RDS itself is written by pyreadr (librdata), so no R installation is
    # needed. Compression is done here rather than by pyreadr: pyreadr's gzip
    # header carries the current mtime, so the same data would give a different
    # file on every run and the nf-test md5 snapshots could never match. Writing
    # the frame uncompressed and recompressing with mtime=0 keeps the output
    # byte-identical across runs. Level 6 is the knee of the curve: level 9 is
    # ~4x slower for ~3% smaller output.
    uncompressed_rds = os.path.join(output_dir, base + ".tmp.rds")
    try:
        pyreadr.write_rds(uncompressed_rds, data)
        with open(uncompressed_rds, "rb") as raw, \
                gzip.GzipFile(filename=args.output, mode="wb", compresslevel=6,
                              mtime=0) as gz:
            shutil.copyfileobj(raw, gz, length=1024 * 1024)
    finally:
        if os.path.exists(uncompressed_rds):
            os.remove(uncompressed_rds)

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
