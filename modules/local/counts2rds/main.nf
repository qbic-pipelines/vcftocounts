process COUNTS2RDS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/79/796ae8263fbd0df33c24322acbe985e1851045461d7a8e335ee0e5df4dcc96e6/data'
        : 'community.wave.seqera.io/library/pandas_rpy2:50401fa0182d5917'}"

    input:
    tuple val(meta), path(csvs)
    tuple val(meta2), path(sample_metadata), path(population_info)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    tuple val(meta), path("*.csv"), emit: csv, optional: true
    tuple val(meta), path("*_dropped_samples.csv"), emit: dropped
    tuple val(meta), path("*_summary.csv"), emit: summary
    tuple val("${task.process}"), val('counts2rds'), eval("echo 1.0.0"), topic: versions, emit: versions_counts2rds

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def write_csv = "${meta.write_csv}" ? "--write-csv" : ""
    """
    counts2rds.py \\
        --input-dir '.' \\
        ${write_csv} \\
        --metadata ${sample_metadata} \\
        --population ${population_info} \\
        --output ${prefix}.rds
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def write_csv = "${meta.csv}" ? "touch ${prefix}.csv" : ""
    """
    touch ${prefix}.rds
    touch ${prefix}_dropped_samples.csv
    touch ${prefix}_summary.csv
    ${write_csv}
    """
}
