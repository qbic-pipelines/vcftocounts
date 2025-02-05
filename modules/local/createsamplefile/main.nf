process CREATE_SAMPLE_FILE {
    tag "$meta.id"
    label 'process_single'

    input:
    val(meta)

    output:
    tuple val(meta), path("*.txt"), emit: samplefile
    path  "versions.yml"          , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    echo "${prefix}" > ${prefix}.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_sample_file: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.0.0'
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_sample_file: $VERSION
    END_VERSIONS
    """
}