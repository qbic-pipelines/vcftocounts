process DIVANALYSIS_BATCHGEN {
    tag "${meta.id}"

    // Generation is single-threaded by design (see hill-batchgen.R), so it takes
    // one CPU. Memory is the binding constraint: it holds the full count table
    // and builds the per-sample-set variant samples. process_medium is the label
    // the old monolith used for the same work, and remains right.
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6c/6c2dd8fc4240adf343ad71f9c56158d87f28b2577f2a6e114b7ab8406f0c4672/data'
        : 'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-future_pruned:2ea4e0e6f556ea1b'}"

    input:
    tuple val(meta), path(rds)

    output:
    // The batch files are the fan-out unit: one entry per batch, so the scheduler
    // can run each batch (and each sample set) as its own task.
    tuple val(meta), path("chunks/chunk_downsample_*_batch_*.rds"), emit: batches
    // The index is what the aggregate step reads to order replicates. It must be
    // published alongside the batches because it carries samples_to_keep.
    tuple val(meta), path("chunks/chunk_*.RData"), emit: chunk_index
    path("manifest.rds"), emit: manifest
    tuple val("${task.process}"), val('batchgen'), eval("echo 1.0.0"), topic: versions, emit: versions_batchgen

    when:
    task.ext.when == null || task.ext.when

    script:
    // Population and superpopulation splits are on by default, matching the hill
    // module, so the per-population levels exist for downstream comparison.
    def evaluatePopulation      = (meta.pop     != null && !meta.pop)     ? '' : '--pop'
    def evaluateSuperPopulation = (meta.superpop != null && !meta.superpop) ? '' : '--superpop'
    def sampleSets   = meta.sample_sets  ?: 10
    def variantSets  = meta.variant_sets ?: 10
    def seed         = meta.seed         ?: 1234
    // hill-lib.R is staged here so the module is self-contained; the scripts find
    // it next to themselves.
    """
    hill-batchgen.R \\
        ${rds} \\
        --chunks chunks \\
        --output manifest.rds \\
        --sample-sets ${sampleSets} \\
        --variant-sets ${variantSets} \\
        --q ${meta.q ?: 1} \\
        --seed ${seed} \\
        ${evaluatePopulation} \\
        ${evaluateSuperPopulation} \\
        ${task.ext.args ?: ''}
    """

    stub:
    """
    mkdir -p chunks
    touch chunks/chunk_1.RData chunks/chunk_downsample_1_batch_1.rds manifest.rds
    """
}
