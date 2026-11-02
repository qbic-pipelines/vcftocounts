/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                               } from '../modules/nf-core/multiqc'
include { TABIX_TABIX                           } from '../modules/nf-core/tabix/tabix'
include { GATK4_GENOTYPEGVCFS                   } from '../modules/nf-core/gatk4/genotypegvcfs'
include { BCFTOOLS_CONCAT                       } from '../modules/nf-core/bcftools/concat'
include { CREATE_SAMPLE_FILE                    } from '../modules/local/createsamplefile'
include { BCFTOOLS_REHEADER                     } from '../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_STATS as BEFORE_FILTER_STATS } from '../modules/nf-core/bcftools/stats'
include { BCFTOOLS_STATS as AFTER_FILTER_STATS  } from '../modules/nf-core/bcftools/stats'
include { BCFTOOLS_VIEW                         } from '../modules/nf-core/bcftools/view'
include { BCFTOOLS_MERGE                        } from '../modules/nf-core/bcftools/merge'
include { BCFTOOLS_ANNOTATE                     } from '../modules/nf-core/bcftools/annotate'
include { VCF2COUNTS                            } from '../modules/local/vcf2counts'
include { COUNTS2RDS                            } from '../modules/local/counts2rds'
include { DIVANALYSIS_BATCHGEN                  } from '../modules/local/divanalysis/batchgen'
include { DIVANALYSIS_BATCHPROC                 } from '../modules/local/divanalysis/batchproc'
include { DIVANALYSIS_AGGREGATE                 } from '../modules/local/divanalysis/aggregate'
// DIVANALYSIS_HILL is kept on disk (the monolith hill-fixed.R still works and is
// the reference the split is validated against) but is no longer wired in: the
// batchgen/batchproc/aggregate path supersedes it. Re-add this include and call
// it directly to fall back to the single-task behaviour.
// include { DIVANALYSIS_HILL                      } from '../modules/local/divanalysis/hill'
include { DIVANALYSIS_PIPECOMP                  } from '../modules/local/divanalysis/pipecomp'
include { DIVANALYSIS_PLOTS                     } from '../modules/local/divanalysis/plots'
include { paramsSummaryMap                      } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_vcftocounts_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VCFTOCOUNTS {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta
    fai
    dict
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:
    ch_multiqc_files = channel.empty()

    //
    // add index to non-indexed VCFs
    //
    (ch_has_index, ch_has_no_index) = ch_samplesheet
        .map { meta, vcf, tbi ->
            // Strip common VCF/GVCF and compression extensions reliably
            def cleanedName = vcf.name.replaceFirst(/(\.g\.vcf\.gz|\.gvcf\.gz|\.vcf\.gz|\.g\.vcf|\.gvcf|\.vcf|\.gz)$/, '')
            [meta + [name: cleanedName], vcf, tbi]
        }
        .branch { _meta, _vcf, tbi ->
            has_index: tbi
            to_index: !tbi
        }

    TABIX_TABIX(ch_has_no_index.map { meta, vcf, _tbi -> [meta, vcf, [], []] })

    ch_indexed = ch_has_no_index
        .join(TABIX_TABIX.out.index)
        .map { meta, vcf, _placeholder, tbi -> [meta, vcf, tbi] }

    // Join both channels back together
    ch_vcf_tbi = ch_has_index.mix(ch_indexed)

    //
    // Convert gvcfs to vcfs
    //
    (ch_gvcf, ch_normal_vcf) = ch_vcf_tbi.branch { it ->
        gvcf: it[0].gvcf
        vcf: !it[0].gvcf
    }


    GATK4_GENOTYPEGVCFS(
        ch_gvcf.map { it -> [it[0], it[1], it[2], [], []] },
        fasta.map { it -> [[id: it.baseName], it] },
        fai.map { it -> [[id: it.baseName], it] },
        dict.map { it -> [[id: it.baseName], it] },
        [[], []],
        [[], []],
    )

    ch_vcf_index = GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi)

    ch_vcf = ch_normal_vcf.mix(ch_vcf_index)

    if (params.filter != null) {

        BEFORE_FILTER_STATS(
            ch_vcf,
            [[], []],
            [[], []],
            [[], []],
            [[], []],
            [[], []],
        )

        ch_vcf_counted = ch_vcf
            .join(BEFORE_FILTER_STATS.out.stats)
            .map { meta, vcf, tbi, stats_file ->
                // Read the file as a string
                def statsText = stats_file.text
                // Use a regex to find the line and extract the number
                def matcher = statsText =~ /number of records:\s*(\d+)/
                def numVar = matcher.find() ? matcher.group(1) : null
                [meta + [numVarBefore: numVar], vcf, tbi]
            }

        //
        // Filter VCFs for pattern given in params.filter
        // Use meta.name for the VCF filename
        //
        BCFTOOLS_VIEW(
            ch_vcf_counted,
            [],
            [],
            [],
        )

        ch_filtered_vcf = BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi)

        AFTER_FILTER_STATS(
            ch_filtered_vcf,
            [[], []],
            [[], []],
            [[], []],
            [[], []],
            [[], []],
        )

        ch_filtered_vcf_counted = ch_filtered_vcf
            .join(AFTER_FILTER_STATS.out.stats)
            .map { meta, vcf, tbi, stats_file ->
                // Read the file as a string
                def statsText = stats_file.text
                // Use a regex to find the line and extract the number
                def matcher = statsText =~ /number of records:\s*(\d+)/
                def numVar = matcher.find() ? matcher.group(1) : null
                [meta + [numVarAfter: numVar, status: "filtered"], vcf, tbi]
            }

        ch_filtered_vcf_counted.collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/bcftools/view/csv") { meta, _vcf, _tbi ->
            def sample = meta.id
            def label = meta.label
            def numVarBefore = meta.numVarBefore
            def numVarAfter = meta.numVarAfter
            def fracFiltered = numVarBefore ? (numVarAfter.toFloat() / numVarBefore.toFloat()) : 0
            def vcf_path = "bcftools/view/${meta.label}/${meta.name}.filter.vcf.gz"
            def tbi_path = "bcftools/view/${meta.label}/${meta.name}.filter.vcf.gz.tbi"

            ["info_filtered_variants.csv", "sample,label,numVarBefore,numVarAfter,fractionFilteredVar,vcf,tbi\n${sample},${label},${numVarBefore},${numVarAfter},${fracFiltered},${vcf_path},${tbi_path}\n"]
        }
    }

    // Add status info to non-filtered VCFs
    ch_vcf = ch_vcf.map {
        meta, vcf, tbi ->
        [meta + [status: "full"], vcf, tbi]
    }
    ch_vcf_all = ch_vcf.mix(ch_filtered_vcf_counted)

    //
    // Concatenate converted VCFs if the entries for "id" "label" and "status" are the same
    //
    (ch_single_vcf, ch_multiple_vcf) = ch_vcf_all
        .map { meta, vcf, tbi ->
            [[meta.id, meta.label, meta.status], meta, vcf, tbi]
        }
        .groupTuple(by: 0)
        .map { _id, metas, vcfs, tbis ->
            def meta = metas[0].findAll { k, _v -> k != 'name' }
            // Take the first meta without filename, they should all be the same for a given ID
            [meta, vcfs.flatten(), tbis.flatten()]
        }
        .branch { it ->
            single: it[1].size == 1
            multiple: it[1].size > 1
        }

    BCFTOOLS_CONCAT(ch_multiple_vcf)

    ch_vcf_index = BCFTOOLS_CONCAT.out.vcf.join(BCFTOOLS_CONCAT.out.tbi)

    ch_vcf_concat = ch_single_vcf.mix(ch_vcf_index)

    if (params.rename) {
        // Create the sample file and add them to ch_vcf_concat
        CREATE_SAMPLE_FILE(ch_vcf_concat.map { it -> it[0] })

        ch_vcf_sample = ch_vcf_concat
            .join(CREATE_SAMPLE_FILE.out.samplefile)
            .map { meta, vcf, _tbi, samplefile -> [meta, vcf, [], samplefile] }

        // Rename samples in vcf with the label
        BCFTOOLS_REHEADER(
            ch_vcf_sample,
            [[], []],
        )

        ch_vcf_index_rh = BCFTOOLS_REHEADER.out.vcf.join(BCFTOOLS_REHEADER.out.index)
    }
    else {
        ch_vcf_index_rh = ch_vcf_concat
    }

    //
    // Merge multiple VCFs per sample (patient) with BCFTOOLS_MERGE
    //

    // Bring all vcfs from one sample into a channel
    // Branch based on the number of VCFs per sample
    //
    // Group by [id, status], NOT by id alone. The two arms (filtered / full) are
    // distinct count matrices downstream, so merging them together destroys the
    // comparison: bcftools merge would fold a sample's filtered and full VCFs
    // into one file, the surviving meta would carry whichever status happened to
    // land first, and the later groupTuple(by: status) would then see a single
    // status across every sample -- collapsing the run to one arm (observed:
    // only all_samples_null produced, all_samples_filtered never emitted).
    (ch_single_id, ch_multiple_id) = ch_vcf_index_rh
        .map { meta, vcf, tbi ->
            [meta.id, meta.status, meta, vcf, tbi]
        }
        .groupTuple(by: [0, 1])
        .map { _id, _status, metas, vcfs, tbis ->
            def meta = metas[0].findAll { k, _v -> k != 'name' }
            // Take the first meta without filename, they should all be the same for a given ID
            [meta, vcfs.flatten(), tbis.flatten()]
        }
        .branch { it ->
            single: it[1].size == 1
            multiple: it[1].size > 1
        }

    // Run BCFTOOLS_MERGE only on samples with multiple VCFs
    BCFTOOLS_MERGE(
        ch_multiple_id.map { meta, vcfs, tbis -> [meta, vcfs, tbis, []] },
        [[], [], []],
    )

    // Merge the results back into a single channel
    ch_merged_vcfs = ch_single_id.mix(
        BCFTOOLS_MERGE.out.vcf.join(BCFTOOLS_MERGE.out.index)
    )

    //
    // remove any IDs from the ID column of the VCF
    //
    if (params.removeIDs) {

        BCFTOOLS_ANNOTATE(
            ch_merged_vcfs.map { it -> [it[0], it[1], it[2], [], [], [], [], []] }
        )

        ch_removedIDs_vcfs = BCFTOOLS_ANNOTATE.out.vcf
    }
    else {
        ch_removedIDs_vcfs = ch_merged_vcfs
    }


    //
    // Convert VCFs to Count Matrices
    //
    VCF2COUNTS(ch_removedIDs_vcfs.map { it -> [it[0], it[1]] })
    ch_metadata_table = channel.fromPath("/mnt/volume/demo-data/balancedsamplesmetadata.csv").collect().map{ file -> [[id:'metadata'], file]}

    ch_counts_csvs = VCF2COUNTS.out.csv
        .map { meta, csv ->
            [meta.status, meta, csv]
        }
        .groupTuple(by: 0)
        .map { status, _metas, csvs ->
            def mode = status == 'full' ? 'null' : 'filtered'
            [[id: "all_samples_${mode}", group: 'all_samples', mode: mode], csvs.flatten()]
        }

    COUNTS2RDS(
        ch_counts_csvs,
        ch_metadata_table
    )

    DIVANALYSIS_BATCHGEN(
        COUNTS2RDS.out.rds
    )

    //
    // Fan out over the generated batches. The generator emits the batch files
    // themselves, so Nextflow schedules them independently instead of looping
    // inside one task. This replaces the monolith's "sequential sample sets,
    // parallel variant sets" strategy: the parallelism is now between tasks,
    // where memory is bounded per task and a single failure costs one batch
    // rather than the whole run.
    //
    // The flatten() is load-bearing. path("chunks/*_batch_*.rds") collects every
    // batch into ONE channel item holding a list of 10 files, which would run a
    // single task over all of them -- exactly what this split exists to avoid.
    // Binding each file as its own queue item gives one task per batch.
    //
    ch_batches = DIVANALYSIS_BATCHGEN.out.batches
        .flatMap { meta, batch_files ->
            // batch_files may be a Path or a list of them depending on how many
            // the generator emitted; normalise before flattening.
            (batch_files instanceof List ? batch_files : [batch_files]).collect { bf ->
                [meta, bf]
            }
        }

    DIVANALYSIS_BATCHPROC(
        ch_batches
    )

    //
    // Reduce, ONE aggregate task PER ARM.
    //
    // The two arms are separate count matrices and must not be merged. Aggregating
    // them together would interleave filtered and null replicates into a single
    // table and destroy the comparison pipecomp.R exists to make, so the arm is
    // carried on the channel key all the way here (
    // meta.id is already arm-specific: all_samples_null / all_samples_filtered).
    //
    // The join is deliberate: passing the results and the counts as two separate
    // process inputs would make Nextflow pair them as a cartesian product (four
    // tasks instead of two). One explicitly joined channel gives exactly one task
    // per arm.
    //
    ch_aggregate_input = DIVANALYSIS_BATCHPROC.out.results
        .map { meta, result -> [meta.id, meta, result] }
        .groupTuple(by: 0)
        .map { id, metas, results -> [id, metas[0], results] }
        .join(
            COUNTS2RDS.out.rds.map { meta, rds -> [meta.id, rds] }
        )
        .map { id, meta, results, counts_rds -> [meta, results, counts_rds] }

    DIVANALYSIS_AGGREGATE(
        ch_aggregate_input
    )

    ch_hill_pair = DIVANALYSIS_AGGREGATE.out.rds
        .map { meta, rds -> [meta.group, meta, rds] }        .groupTuple(by: 0)
        .map { _group, metas, rds ->
            def iNull = metas.findIndexOf { m -> m.mode == 'null' }
            def iFiltered = metas.findIndexOf { m -> m.mode == 'filtered' }
            // Without --filter there is no filtered arm, so both sides take the null rds;
            // pipecomp.R warns when --filtered and --null resolve to the same file.
            def nullRds = iNull >= 0 ? rds[iNull] : rds[0]
            def filteredRds = iFiltered >= 0 ? rds[iFiltered] : nullRds
            [[id: 'all_samples'], filteredRds, nullRds]
        }

    // Only getting numbers and displaying them
    DIVANALYSIS_PIPECOMP(
        ch_hill_pair
    )

    // PLOTS needs both arms: null and filtered.  The aggregate already groups
    // by arm, so extract the null and filtered rds the same way as pipecomp does.
    def plots_pair = ch_hill_pair
        .map { meta, filteredRds, nullRds -> [meta, nullRds, filteredRds] }

    DIVANALYSIS_PLOTS(
        plots_pair
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'vcftocounts_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'vcftocounts'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
}
