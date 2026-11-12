# qbic-pipelines/vcftocounts: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Tabix](#tabix) - Recompresses and indexes (g.)vcf files
- [GenotypeGVCFs](#genotypegvcfs) - Converts g.vcf files to vcf with GATK
- [Filter VCFs](#filter-vcfs) - Filters the VCF based on a string given to the `filter` param with bcftools/view
- [Concatenate VCFs](#concatenate-vcfs) - Concatenates all vcfs that have the same id and the same label with bcftools/concat
- [Rename Samples](#rename-samples) - Changes the sample name in the vcf file to the label with bcftools/reheader
- [Merge VCFs](#merge-vcfs) - Merges all vcfs from the same sample with bcftools/merge
- [Remove IDs](#remove-ids) - Removes entries in ID field with bcftools/annotate
- [Convert to matrix](#convert-to-matrix) - Converts the (merged) vcfs to a matrix using a custom Python script written by @ellisdoro
- [Merge to RDS](#merge-to-rds) - Merges the per-sample matrices into one RDS with sample metadata using `counts2rds`
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Tabix

<details markdown="1">
<summary>Output files</summary>

- `htslib/{meta.label}/`
  - `{filename}.indexed.vcf.gz`: bgzf-compressed vcf file.
  - `{filename}.indexed.vcf.gz.tbi`: tabix index of the vcf file.

</details>

For all `(g).vcf` files that are given to the pipeline without an index, `htslib/bgziptabix` recompresses the input to bgzf (if it is not already bgzf-compressed) and generates a tabix index with `.tbi` extension. The indexed VCF and its index are emitted together by the same process, so they always belong to each other. The suffix `.indexed` keeps the output distinct from the input file name.

### GenotypeGVCFs

<details markdown="1">
<summary>Output files</summary>

- `gvcftovcf/{meta.label}/`
  - `{filename}.normal.vcf.gz`: normal vcf file based on gvcf input.
  - `{filename}.normal.vcf.gz.tbi`: tabix index of the vcf file.

</details>

The GATK GenotypeGVCFs module translates genotype (g) vcf files into classic vcf files. The key difference between a regular VCF and a GVCF is that the GVCF has records for all sites, whether there is a variant call there or not.

### Filter VCFs

<details markdown="1">
<summary>Output files</summary>

- `bcftools/view/{meta.label}/`
  - `{filename}.filter.vcf.gz`: vcf file with filtered variants.
  - `{filename}.filter.vcf.gz.tbi`: tabix index of the vcf file.
- `bcftools/view/csv/`
  - `info_filtered_variants.csv`: table containing numbers of variants before and after filtering and the fraction.
- `bcftools/stats/{meta.label}/`
  - `{filename}.bcftools_stats.txt`: bcftools/stats output on vcf before filtering.
  - `{filename}.filter.bcftools_stats.txt`: bcftools/stats output on vcf after filtering.

</details>

VEP annotated VCF files can be filtered for certain flags present after VEP annotation. Notably, this enables filtering for variants with certain impact levels or consequences. Filtering will extract variants that match the specific patterns. When filtering is turned on, `bcftools/stats` is run twice - before and after filtering. The resulting stats are used to compile a CSV file (`info_filtered_variants.csv`) with variant counts and the fraction of variants which are kept.

### Concatenate VCFs

<details markdown="1">
<summary>Output files</summary>

- `bcftools/concat/`
  - `{sample}.{label}.concat.vcf.gz`: vcf file containing all variants from files with same {sample} and {label}.
  - `{sample}.{label}.concat.vcf.gz.tbi`: tabix index of the vcf file.

</details>

Some variant calling pipelines will return multiple (g)VCF files for one patient. The `concatenate` function of `bcftools` is used to add these VCFs to one VCF.

### Rename Samples

<details markdown="1">
<summary>Output files</summary>

- `bcftools/reheader/{sample}/`
  - `{sample}.{label}.reheader.vcf.gz`: vcf file with renamed sample {label}.
  - `{sample}.{label}.reheader.vcf.gz.tbi`: tabix index of the vcf file.

</details>

To make enable the comparison of the finalized CSV files, `bcftools reheader` can be enabled to rename the variant sample name from the generic name given by the variant caller to a custom label given with the samplesheet. Can be turned off with `--rename false`.

### Merge VCFs

<details markdown="1">
<summary>Output files</summary>

- `bcftools/merge/`
  - `{sample}.merge.vcf.gz`: merged vcf file with multiple samples (one per pipeline/caller).
  - `{sample}.merge.vcf.gz.tbi`: tabix index of the vcf file.

</details>

To enable comparison of different variant callers or variant calling pipelines, all VCFs that come from the same sample are merged based on the sample ID submitted by the user.

### Remove IDs

<details markdown="1">
<summary>Output files</summary>

- `bcftools/annotate/`
  - `{sample}.IDremoved.vcf.gz`: vcf file without rsIDs as variant names.
  - `{sample}.IDremoved.vcf.gz.tbi`: tabix index of the vcf file.

</details>

Removes entries in the `ID` column of the VCF using `bcftools annotate -x ID` to prepare for matrix conversion. If the entries are not removed, the matrix conversion script will use available IDs instead of chromosome + position to map the variants. Can be turned off with `--removeIDs false`.

### Convert to matrix

<details markdown="1">
<summary>Output files</summary>

- `vcf2counts/`
  - `{sample}.csv`: csv file containing the variants, one column per {label}

</details>

A custom Python script is used to convert the finalized VCF to a CSV which can be used for further downstream analysis. Script was written by [Dorothy Ellis](https://github.com/ellisdoro). The VCF is streamed line by line rather than loaded into memory, so large cohort VCFs can be converted without a memory penalty.

The CSV holds one row per variant, labelled `{chrom}:{pos}_{REF}/{ALT}` using the first ALT allele only, and one column per sample. Cells encode the genotype: `0` for hom-ref (`0/0`), `1` for het (`0/1`, `1/0`), `2` for hom-alt (`1/1`) and `NA` for anything else, such as multi-allelic genotypes (`0/2`).

The first column holds the variant id (chromosome and position, unless rsIDs were kept), and every further column is one caller, holding the genotype call per variant as `0` (homozygous reference), `1` (heterozygous), `2` (homozygous alternate) or `NA` (not called). `NA` is used for both a genotype the caller left uncalled and a site the caller did not report at all, so `0` always means a real homozygous-reference call. The column name is the `label` from the samplesheet, or the sample name recorded in the VCF when `--rename false` is used. A sample that had no variants surviving filtering still produces a header-only CSV rather than no file at all.

### Merge to RDS

<details markdown="1">
<summary>Output files</summary>

- `counts2rds/`
  - `all_samples.rds`: one R object holding all samples, all callers, and the sample metadata. This is the file the downstream diversity analysis reads.
  - `all_samples.csv`: csv copy of the same table, written only when `write_csv` is set, because the RDS is binary and not inspectable by eye.
  - `all_samples_summary.csv`: run metrics as `metric,value` rows (`n_samples`, `n_callers`, `callers`, `n_rows`, `n_variants`, `n_duplicate_keys_collapsed`, `n_samples_dropped`, `n_rows_without_metadata`, `n_na_genotypes`, `n_populations`, `n_superpopulations`, `metadata_file`, `population_lookup_file`).
  - `all_samples_dropped_samples.csv`: samples that were dropped, with `sample_name`, `source_file` and `reason`. Written even when no sample was dropped, so the file always exists.

</details>

The per-sample CSVs from the previous step are pooled into a single table and saved as an RDS. The table has one row per sample and variant, with these columns:

| Column                         | Description                                      |
| ------------------------------ | ------------------------------------------------ |
| `sample_name`                  | Sample id, taken from the CSV filename.          |
| `variant_called`               | Variant id.                                      |
| one column per `label`         | Genotype call per caller: `0`, `1`, `2` or `NA`. |
| `population_code`              | Population code, e.g. `GBR`.                     |
| `super_population_code`        | Superpopulation code, e.g. `EUR`.                |
| `population_description`       | Population description, e.g. `British`.          |
| `super_population_description` | Superpopulation description, e.g. `European`.    |

The four metadata columns are filled from `--sample_metadata` and `--population_info`; see the [usage documentation](usage.md#sample-metadata) for the accepted file shapes. Three behaviours are worth knowing before reading the RDS:

- The genotype columns hold `0` (homozygous reference), `1` (heterozygous), `2` (homozygous alternate) and `NA` (not called). `NA` is a real missing value, not a stand-in for `0`: it covers both a variant that is absent from that sample's `vcf2counts` CSV, because the caller did not report it there, and a genotype `vcf2counts` itself marked `NA` (uncalled or multi-allelic). Read it with `is.na()`, and pass `na.rm = TRUE` to any arithmetic that should ignore those cells. The count is reported as `n_na_genotypes` in the summary and in the task log.

  A hom-ref `0` and a not-called `NA` answer different questions. `0` is a call: the site was evaluated for that sample and the genotype was reference, so that sample does not carry this particular ALT allele. `NA` claims nothing. A variant row exists in the VCF because at least one sample carries the ALT allele, and in that same row every other sample is `0/0`, so `0` here means "this sample does not have the variant, the other sample does". That is a genotype, not an absence of information.

  `NA` appears where that reasoning breaks down. When the per-sample CSV has no row for a variant at all, the caller never reported the site for that sample and nothing can be inferred; this is the common case in a multi-caller run, because each caller emits its own variant set and the union across samples is larger than any single one. It also appears for a genotype `vcf2counts` marked `NA`, i.e. uncalled (`./.`) or multi-allelic (`0/2`), where the call exists but is not one of the three states the matrix can express. Filling either with `0` would assert a homozygous-reference genotype nobody reported.

- Repeated `(sample_name, variant_called)` pairs are collapsed by taking the maximum per caller column, so the strongest genotype call wins. `data.table` rejects duplicate row names and the downstream code builds rownames as `sample_name` + `_` + `variant_called`.
- Rows whose sample has no metadata match are kept, but their metadata columns are `NA` and they drop out of any per-population downstream stage. The count is reported as `n_rows_without_metadata` in the summary and as a warning in the task log.

The RDS is written by [`pyreadr`](https://ofajardo.github.io/pyreadr/), so the module needs no R installation. Count columns are stored as R integers and the metadata columns as character with `NA`, and the file is gzip-compressed without a timestamp, so the same input always gives a byte-identical RDS.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
