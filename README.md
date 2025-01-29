# qbic-pipelines/vcftomat

[![GitHub Actions CI Status](https://github.com/qbic-pipelines/vcftomat/actions/workflows/ci.yml/badge.svg)](https://github.com/qbic-pipelines/vcftomat/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/qbic-pipelines/vcftomat/actions/workflows/linting.yml/badge.svg)](https://github.com/qbic-pipelines/vcftomat/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.14616650-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.14616650)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/qbic-pipelines/vcftomat)

## Introduction

**qbic-pipelines/vcftomat** is a bioinformatics pipeline that ...

1. Indexes (g.)vcf files ([`tabix`](http://www.htslib.org/doc/tabix.html))
2. Converts g.vcf files to vcf with `genotypegvcf` ([`GATK`](https://gatk.broadinstitute.org/hc/en-us))
3. Concatenates all vcfs that have the same id and the same label with `bcftools/concat` ([`bcftools`](https://samtools.github.io/bcftools/bcftools.html))
4. Changes the sample name in the vcf file to the filename with `bcftools/reheader` ([`bcftools`](https://samtools.github.io/bcftools/bcftools.html)) - This can be turned off by adding `--rename false` to the `nextflow run` command.
5. Merges all vcfs from the same sample with `bcftools/merge` ([`bcftools`](https://samtools.github.io/bcftools/bcftools.html))
6. Converts the (merged) vcfs to a matrix using a custom R script written by @ellisdoro ([`R`](https://www.r-project.org/))
7. Collects all reports into a MultiQC report ([`MultiQC`](http://multiqc.info/))

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,label,gvcf,vcf_path,vcf_index_path
SAMPLE-1,pipelineA-callerA,false,path/to/vcf.gz,path/to/.vcf.gz.tbi
SAMPLE-1,pipelineB-callerA,false,path/to/vcf.gz,path/to/.vcf.gz.tbi
SAMPLE-2,pipelineB-callerB,true,path/to/g.vcf.gz,path/to/g.vcf.gz.tbi
SAMPLE-2,pipelineB-callerB,true,path/to/g.vcf.gz,path/to/g.vcf.gz.tbi
```

Each row represents a VCF file coming from a sample. The `label` column enables concatenation of vcfs (for example when the pipeline produces different vcfs for chrM and chrY). The `gvcf` column indicates whether the file is a g.vcf file or not. The `vcf_path` and `vcf_index_path` columns contain the path to the VCF file and its index, respectively.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run qbic-pipelines/vcftomat \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

qbic-pipelines/vcftomat was originally written by Famke Bäuerle, Dorothy Ellis.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use qbic-pipelines/vcftomat for your analysis, please cite it using the following doi: [10.5281/zenodo.14616650](https://doi.org/10.5281/zenodo.14616650)

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
