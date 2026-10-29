process DIVANALYSIS_PLOTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    // Built by Wave from this module's environment.yml, which adds ggplot2 and pals
    // to the shared divanalysis base. This module CANNOT reuse the hill container:
    // that one has neither package, so plots.R would die on library(ggplot2).
    //
    // A single image reference is used deliberately rather than the usual
    // docker/singularity ternary: no singularity blob was emitted for this build, and
    // reusing the hill blob here would point at an image without ggplot2 - a silent
    // failure. With `-profile wave`, Wave resolves this docker reference to the
    // engine in use. Regenerate the image with the wave profile if environment.yml
    // changes.
    container 'community.wave.seqera.io/library/r-base_r-data.table_r-dplyr_r-docopt_pruned:74bdd4f4331138b0'

    input:
    tuple val(meta), path(null_rds), path(filtered_rds)

    output:
    tuple val(meta), path("*_mqc.png"), emit: plots
    tuple val(meta), path("*_nulltest_*.csv"), emit: csv, optional: true
    tuple val(meta), path("*_plots_metadata.rds"), emit: rds
    tuple val("${task.process}"), val('hillplots'), eval("echo 1.0.0"), topic: versions, emit: versions_plots

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Default prefix is the bare id, matching the pipelinediff module: the script adds
    // its own '_plots_' infix to the metadata file, so '${meta.id}_plots' here would
    // produce '<id>_plots_plots_metadata.rds'.
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Splits are on by default, matching the hill module, so the per-population and
    // per-superpopulation densities and tests are produced without extra wiring.
    // Set meta.pop / meta.superpop to false to opt out per run. When the rds has no
    // usable split for a level the script skips just that level.
    def evaluatePopulation = (meta.pop != null && !meta.pop) ? '' : '--pop'
    def evaluateSuperPopulation = (meta.superpop != null && !meta.superpop) ? '' : '--superpop'
    def nPermutations = meta.n_permutations ? "--n-permutations ${meta.n_permutations}" : ''
    def plotSeed = meta.seed ? "--seed ${meta.seed}" : ''
    def qHillNumber = meta.q ? "--q ${meta.q}" : ''
    """
    plots.R \\
        --null ${null_rds} \\
        --filtered ${filtered_rds} \\
        --output ${prefix} \\
        --outdir . \\
        ${evaluatePopulation} \\
        ${evaluateSuperPopulation} \\
        ${nPermutations} \\
        ${plotSeed} \\
        ${qHillNumber} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_001_Null_Overall_TD_beta_mqc.png
    touch ${prefix}_plots_metadata.rds
    """
}
