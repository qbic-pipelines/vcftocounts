process VCF2MAT {
    tag "$meta.id"
    label 'process_single'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://famkebaeuerle/vcf2mat:1.0.0' :
        'famkebaeuerle/vcf2mat:1.0.0' }"

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
