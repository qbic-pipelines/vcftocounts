process COUNTS2RDS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6a/6a5f40d32ae6ef6abebdb9f2898f84f31f52fe8bee37b71f9b8ff97b1ebf1d35/data'
        : 'community.wave.seqera.io/library/r-base_r-data.table_r-optparse:bdc0c95ecea704b8'}"

    input:
    tuple val(meta), path(csvs, stageAs: 'counts/*')
    tuple val(meta2), path(sample_metadata), path(population_info)

    output:
    tuple val(meta), path("${prefix}.null.rds"),            emit: nulltable
    tuple val(meta), path("${prefix}.variants.rds"),        emit: variants
    tuple val(meta), path("${prefix}.info.rds"),            emit: info
    // Optional on the Nextflow side only because of the call-pattern ceiling below.
    tuple val(meta), path("${prefix}.q2.rds"),              emit: q2,      optional: true
    tuple val(meta), path("${prefix}.summary.csv"),         emit: summary
    // Kept identical to the other local modules (emit before topic): this pattern already runs
    // in vcf2counts and createsamplefile.
    tuple val("${task.process}"), val('counts2rds'), eval("echo 1.0.0"), topic: versions, emit: versions_counts2rds

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def args       = task.ext.args ?: ''
    def population = population_info ? "--population ${population_info}" : ''
    """
    counts2rds.R \\
        --counts counts/ \\
        --metadata ${sample_metadata} \\
        ${population} \\
        --threads ${task.cpus ?: 1} \\
        ${args} \\
        --out-dir . \\
        --prefix ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.null.rds ${prefix}.variants.rds ${prefix}.info.rds
    touch ${prefix}.q2.rds ${prefix}.summary.csv
    """
}
