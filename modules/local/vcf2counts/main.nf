process VCF2COUNTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2e/2ee25e965bee9d9deaec844671537eee86cc71fd6fa2f9513482bd45b2235921/data'
        : 'community.wave.seqera.io/library/python:3.14.8--7620d6d80d3ba0cb'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    tuple val("${task.process}"), val('vcf2counts'), eval("echo 2.0.0"), topic: versions, emit: versions_vcf2counts

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf2counts.py \\
        --output ${prefix}.csv \\
        ${vcf}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    """
}
