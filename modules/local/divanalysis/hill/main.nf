process DIVANALYSIS_HILL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6c/6c2dd8fc4240adf343ad71f9c56158d87f28b2577f2a6e114b7ab8406f0c4672/data'
        : 'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-future_pruned:2ea4e0e6f556ea1b'}"

    input:
    tuple val(meta), path(rds)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    tuple val("${task.process}"), val('hillcalculation'), eval("echo 1.0.0"), topic: versions, emit: versions_hill

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_hill"
    // Population and superpopulation splits are produced by default so the
    // per-population/per-superpopulation levels are available downstream without
    // every caller having to remember the flags. main.R leaves these off; here
    // they are on because the split is the point of the downstream comparison.
    // Enabling them does not change the total result: the flags only create the
    // split indices alongside, they draw no extra randomness, so the total
    // PairwiseTotalData is bit-identical either way (verified).
    // meta.pop / meta.superpop can be set to false per run to opt out.
    def evaluatePopulation      = (meta.pop     != null && !meta.pop)     ? '' : '--pop'
    def evaluateSuperPopulation = (meta.superpop != null && !meta.superpop) ? '' : '--superpop'
    """
    hill-fixed.R \\
        --mode ${meta.mode} \\
        --cores ${task.cpus} \\
        ${evaluatePopulation} \\
        ${evaluateSuperPopulation} \\
        ${task.ext.args ?: ''} \\
        --output ${prefix}.rds \\
        ${rds}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rds
    """
}
