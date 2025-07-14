# qbic-pipelines/vcftocounts: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.2dev

### `Added`

- [#34](https://github.com/qbic-pipelines/vcftocounts/pull/34) - Swap CI tests to nf-test and fix small channel issue
- [#40](https://github.com/qbic-pipelines/vcftocounts/pull/40) - Create csv containing info on filtered variants

### `Fixed`

- [#33](https://github.com/qbic-pipelines/vcftocounts/pull/33) - Back to dev 2.0.2dev
- [#35](https://github.com/qbic-pipelines/vcftocounts/pull/35) - Template update to version 3.2.1
- [#36](https://github.com/qbic-pipelines/vcftocounts/pull/36) - Update all nf-core modules
- [#37](https://github.com/qbic-pipelines/vcftocounts/pull/37) - Fix output of filtering module to avoid overwriting
- [#38](https://github.com/qbic-pipelines/vcftocounts/pull/38) - Template update to version 3.3.1 + update multiqc
- [#41](https://github.com/qbic-pipelines/vcftocounts/pull/41) - Template update to version 3.3.2

### `Dependencies`

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| BCFtools   | 1.2         | 1.21        |
| Tabix      | 1.2         | 1.21        |
| MultiQC    | 1.27        | 1.29        |

## v2.0.1 - Pipe Cerulean - 16.04.2025

### `Fixed`

- [#29](https://github.com/qbic-pipelines/vcftocounts/pull/29) - Use whole name instead of baseName for gvcf to vcf
- [#30](https://github.com/qbic-pipelines/vcftocounts/pull/30) - Use filtered vcfs for subsequent processes
- [#31](https://github.com/qbic-pipelines/vcftocounts/pull/31) - Prepare release 2.0.1

## v2.0.0 - Rad Sepia - 12.02.2025

Initial release of (renamed) qbic-pipelines/vcftocounts.

### `Changed`

- [#25](https://github.com/qbic-pipelines/vcftocounts/pull/25) - Rename pipeline to vcftocounts + Prepare Release 2.0.0
- [#28](https://github.com/qbic-pipelines/vcftocounts/pull/28) - Move test data to test-datasets, create filter nf-test

### `Added`

- [#22](https://github.com/qbic-pipelines/vcftocounts/pull/22) - Remove ID column entries from VCFs
- [#23](https://github.com/qbic-pipelines/vcftocounts/pull/23) - Add filtering of VEP annotated VCF files using bcftools/view (no conda for NXF_VER <= 24.10.2)

### `Fixed`

- [#15](https://github.com/qbic-pipelines/vcftocounts/pull/15) - template update to v3.1.1
- [#16](https://github.com/qbic-pipelines/vcftocounts/pull/16) - Fix bcftools/reheader bug
- [#21](https://github.com/qbic-pipelines/vcftocounts/pull/21) - template update to v3.2.0
- [#24](https://github.com/qbic-pipelines/vcftocounts/pull/24) - Change branching logic to enable caching again + rename vcf2mat to vcf2counts (thank you @awgymer for helping)

## v1.1.0 - Newton Puccoon - 08.01.2025

### `Added`

- [#7](https://github.com/qbic-pipelines/vcftocounts/pull/7) - samplenames to columns
- [#8](https://github.com/qbic-pipelines/vcftocounts/pull/8) - concat for sample, label pairs

### `Fixed`

- [#5](https://github.com/qbic-pipelines/vcftocounts/pull/5) - filename collision
- [#10](https://github.com/qbic-pipelines/vcftocounts/pull/10) - prepare release 1.1.0

## v1.0.0 - Curie Purpureal - 16.12.2024

Initial release of qbic-pipelines/vcftocounts, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
