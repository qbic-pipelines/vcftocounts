process DIVANALYSIS_AGGREGATE {
    tag "${meta.id}"

    // The reduce step: cheap next to generation and processing, but it holds the
    // merged replicate table in memory, so it is not process_single. The original
    // did this work inside the monolith's task; here it is its own stage.
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6c/6c2dd8fc4240adf343ad71f9c56158d87f28b2577f2a6e114b7ab8406f0c4672/data'
        : 'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-future_pruned:2ea4e0e6f556ea1b'}"

    // Aggregation is all-or-nothing: a partial table silently changes every
    // downstream comparison, so retry hard and fail loudly rather than publish
    // something misleading.
    errorStrategy 'retry'
    maxRetries 3

    input:
    // ONE tuple with three elements, not two separate inputs. Two inputs would be
    // paired as a cartesian product by Nextflow (2 arms x 2 counts = 4 tasks);
    // joining upstream keeps it to exactly one task per arm.
    //
    // result_files is every batch result for THIS arm. The upstream grouping
    // means all replicates arrive together, which list_transpose and the
    // reorganise step require -- and why duplicates/holes are a hard error there.
    tuple val(meta), path(result_files, stageAs: 'results/*'), path(rds)

    output:
    // Literal name, not a script-local: output paths resolve before `script:`
    // runs. task.ext.prefix overrides it when a caller wants a different name.
    tuple val(meta), path("*_hill.rds"), emit: rds
    tuple val("${task.process}"), val('hillaggregate'), eval("echo 1.0.0"), topic: versions, emit: versions_aggregate

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_hill"
    def variantSets = meta.variant_sets ?: 10
    // The batchproc stage names each output after its batch file, so the staged
    // path is results/<batch-stem>.rds. Build the argument list explicitly and
    // sort it, so the command line is reproducible and Nextflow's task hash does
    // not depend on filesystem enumeration order.
    def inputs = result_files.sort { it.name }.collect { "results/${it.baseName}.rds" }.join(' ')
    """
    hill-aggregate.R \\
        ${inputs} \\
        --counts ${rds} \\
        --output ${prefix}.rds \\
        --variant-sets ${variantSets} \\
        ${task.ext.args ?: ''}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_hill"
    """
    touch ${prefix}.rds
    """
}
