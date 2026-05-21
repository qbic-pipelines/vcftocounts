process VCF2COUNTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6c/6c2dd8fc4240adf343ad71f9c56158d87f28b2577f2a6e114b7ab8406f0c4672/data'
        : 'community.wave.seqera.io/library/bioconductor-variantannotation_r-docopt_r-matrix:3cf2f20fdc477746'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    tuple val("${task.process}"), val('vcf2counts'), eval("echo 1.0.0"), topic: versions, emit: versions_vcf2counts

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf2counts.R \\
        --output ${prefix}.csv \\
        ${vcf}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    """
}
