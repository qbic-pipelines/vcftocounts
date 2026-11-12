process COUNTS2RDS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9d/9d2d32679d6075e6f6f95600e96531beb321a54e510a1e9ccc472b2174dd0662/data'
        : 'community.wave.seqera.io/library/pandas_pyreadr:65b993d44598cac9'}"

    input:
    tuple val(meta), path(csvs, stageAs: 'counts/*')
    tuple val(meta2), path(sample_metadata), path(population_info)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    tuple val(meta), path("*.csv"), emit: csv, optional: true
    tuple val(meta), path("*_dropped_samples.csv"), emit: dropped
    tuple val(meta), path("*_summary.csv"), emit: summary
    tuple val("${task.process}"), val('counts2rds'), eval("echo 1.0.0"), topic: versions, emit: versions_counts2rds

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def write_csv = meta.write_csv ? "--write-csv" : ""
    // Both metadata files are optional: an unset parameter reaches argparse as
    // the literal string "null", which it rejects.
    def metadata = sample_metadata ? "--metadata ${sample_metadata}" : ""
    def population = population_info ? "--population ${population_info}" : ""
    """
    counts2rds.py \\
        --input-dir counts/ \\
        ${write_csv} \\
        ${metadata} \\
        ${population} \\
        --output ${prefix}.rds
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def write_csv = meta.write_csv ? "touch ${prefix}.csv" : ""
    """
    touch ${prefix}.rds
    touch ${prefix}_dropped_samples.csv
    touch ${prefix}_summary.csv
    ${write_csv}
    """
}
