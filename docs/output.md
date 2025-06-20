# qbic-pipelines/vcftocounts: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Tabix](#tabix) - Indexes (g.)vcf files
- [GenotypeGVCFs](#genotypegvcfs) - Converts g.vcf files to vcf with GATK
- [Filter VCFs](#filter-vcfs) - Filters the VCF based on a string given to the `filter` param with bcftools/view
- [Concatenate VCFs](#concatenate-vcfs) - Concatenates all vcfs that have the same id and the same label with bcftools/concat
- [Rename Samples](#rename-samples) - Changes the sample name in the vcf file to the label with bcftools/reheader
- [Merge VCFs](#merge-vcfs) - Merges all vcfs from the same sample with bcftools/merge
- [Remove IDs](#remove-ids) - Removes entries in ID field with bcftools/annotate
- [Convert to matrix](#convert-to-matrix) - Converts the (merged) vcfs to a matrix using a custom R script written by @ellisdoro
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Tabix

<details markdown="1">
<summary>Output files</summary>

- `tabix`
  - `{filename}.vcf.gz.tbi`: tabix index of the vcf file.

</details>

Tabix generated index files with `.tbi` extension for all `(g).vcf` files that are given to the pipeline without index.

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

VEP annotated VCF files can be filtered for certain flags present after VEP annotation. Notably, this enables filtering for variants with certain impact levels or consequences. Filtering will produces VCF files holding just the variants matching the specific patterns. Running the filtering part also enables `bcftools/stats` which is used to compute the fraction of variants which are kept. These results are found in a csv file called `info_filtered_variants.csv`.

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

Removes entries in the `ID` column of the VCF using `bcftools annotate -x ID` to prepare for matrix conversion. If the entries are not removed, the R script will use available IDs instead of chromosome + position to map the variants. Can be turned off with `--removeIDs false`.

### Convert to matrix

<details markdown="1">
<summary>Output files</summary>

- `vcf2counts/`
  - `{sample}.csv`: csv file containing the variants, one column per {label}

</details>

A custom R script is used to convert the finalized VCF to a CSV which can be used for further downstream analysis. Script was written by [Dorothy Ellis](https://github.com/ellisdoro).

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

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
