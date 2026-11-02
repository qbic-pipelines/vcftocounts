process DIVANALYSIS_BATCHPROC {
    tag "${meta.id}"

    // THE unit of parallel work. Sized for ONE batch (~1.2 GB peak measured),
    // not for a whole sample set, which is the point of the split: several of
    // these fit on a small host at once, where the old monolith needed one task
    // holding everything.
    //
    // process_low rather than process_medium: 12 GB is generous headroom over the
    // measured 1.2 GB, and the batching is what lets this step use a smaller
    // label than the 36 GB the monolith asked for.
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6c/6c2dd8fc4240adf343ad71f9c56158d87f28b2577f2a6e114b7ab8406f0c4672/data'
        : 'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-future_pruned:2ea4e0e6f556ea1b'}"

    // A batch that dies is re-run alone, which is the retry granularity this
    // split buys. Static rather than a dynamic closure: a dynamic errorStrategy
    // has been seen to fail with "Directive `errorStrategy` doesn't support
    // dynamic value (or context not yet initialized)" precisely when the process
    // errors, which is the worst possible time. After the retries are exhausted
    // the run aborts, by design -- the aggregate step refuses holes by default,
    // so a silently partial table is not a possibility.
    errorStrategy 'retry'
    maxRetries 2

    input:
    // meta is carried through for provenance; batch_file is a SINGLE file, which
    // is what makes Nextflow run one task per batch. The workflow flattens the
    // generator's file list to get that (see vcftocounts.nf).
    tuple val(meta), path(batch_file)

    output:
    // Glob rather than interpolate a script-local: output paths are resolved
    // before the script block runs, so `prefix` is not in scope here. The
    // aggregate stage reads whatever names this publishes.
    //
    // The result suffix must NOT reproduce the input's name. Naming the output
    // after the batch file (chunk_downsample_N_batch_1.rds) was doubly broken:
    // (a) Nextflow excludes input file names from output glob matching, so the
    // task always reported "Missing output file(s)" and retried forever, and
    // (b) the input is staged as a SYMLINK, so writing to that name wrote
    // through it and destroyed the generator's batch file in the work dir.
    tuple val(meta), path("*_batch_*.result.rds"), emit: results
    tuple val("${task.process}"), val('batchproc'), eval("echo 1.0.0"), topic: versions, emit: versions_batchproc

    when:
    task.ext.when == null || task.ext.when

    script:
    // Result name = batch stem + ".result", so it can never collide with the
    // input name (see the note on the output block). The batch file itself is
    // never written to, which is what keeps the generator's output intact.
    def prefix = "${batch_file.baseName}.result"
    def variantSets = meta.variant_sets ?: 10
    """
    hill-batchproc.R \\
        --batch ${batch_file} \\
        --output ${prefix}.rds \\
        --variant-sets ${variantSets} \\
        --q ${meta.q ?: 1} \\
        --strict \\
        ${task.ext.args ?: ''}
    """

    stub:
    def prefix = "${batch_file.baseName}.result"
    """
    touch ${prefix}.rds
    """
}
