process CREATE_SAMPLE_FILE {
    tag "${meta.id}"
    label 'process_single'

    input:
    val meta

    output:
    tuple val(meta), path("*.txt"), emit: samplefile
    tuple val("${task.process}"), val('create_sample_file'), eval("echo 1.0.0"), topic: versions, emit: versions_create_sample_file                     

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${prefix}" > ${prefix}.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    """
}
