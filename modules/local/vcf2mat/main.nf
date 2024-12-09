process VCF2MAT {
    tag "$meta.id"
    label 'process_single'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/famkebaeuerle/vcf2mat:1.0.0' :
        'community.wave.seqera.io/library/bioconductor-variantannotation_r-docopt_r-matrix:3cf2f20fdc477746' }"

    input:
    tuple val(meta), path(vcf)


    output:
    tuple val(meta), path("*.csv") , emit: csv
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    vcf2counts.R \\
        --output ${prefix}.csv \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2counts.R: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2counts.R: $VERSION
    END_VERSIONS
    """
}
