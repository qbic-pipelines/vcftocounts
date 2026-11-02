process DIVANALYSIS_PIPECOMP {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6c/6c2dd8fc4240adf343ad71f9c56158d87f28b2577f2a6e114b7ab8406f0c4672/data'
        : 'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-future_pruned:2ea4e0e6f556ea1b'}"

    input:
    tuple val(meta), path(filtered_rds), path(null_rds)

    output:
    tuple val(meta), path("*_filtered_vs_null_*.csv"), emit: csv
    tuple val(meta), path("*_filtered_vs_null_metadata.rds"), emit: rds
    tuple val("${task.process}"), val('pipecomp'), eval("echo 1.0.0"), topic: versions, emit: versions_pipecomp

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    // Defaults ON to match the hill module, so the splits it produces are actually
    // compared. Set meta.pop / meta.superpop to false to opt out per run.
    def evaluatePopulation = (meta.pop != null && !meta.pop) ? '' : '--pop'
    def evaluateSuperPopulation = (meta.superpop != null && !meta.superpop) ? '' : '--superpop'
    def qHillNumber = meta.q ? "--q ${meta.q}" : ''
    def nPermutations = meta.n_permutations ? "--n-permutations ${meta.n_permutations}" : ''
    def seed = meta.seed ? "--seed ${meta.seed}" : ''
    """
    pipecomp.R \\
        --filtered ${filtered_rds} \\
        --null ${null_rds} \\
        --output ${prefix} \\
        --outdir . \\
        ${args} \\
        ${qHillNumber} \\
        ${nPermutations} \\
        ${seed} \\
        ${evaluatePopulation} \\
        ${evaluateSuperPopulation}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_filtered_vs_null_total_pattern_summary.csv
    touch ${prefix}_filtered_vs_null_metadata.rds
    """
}
