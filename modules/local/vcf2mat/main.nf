process VCF2MAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data':
        'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd' }"

    input:
    tuple val(meta), path(vcf)


    output:
    tuple val(meta), path("*.csv"), emit: csv
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    Rscript VCFtoMat.R \\
        --output ${prefix}.csv \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VCFtoMat.R: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VCFtoMat.R: $VERSION
    END_VERSIONS
    """
}
