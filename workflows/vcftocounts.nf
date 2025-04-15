/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { TABIX_TABIX            } from '../modules/nf-core/tabix/tabix/main'
include { GATK4_GENOTYPEGVCFS    } from '../modules/nf-core/gatk4/genotypegvcfs/main'
include { BCFTOOLS_CONCAT        } from '../modules/nf-core/bcftools/concat/main'
include { CREATE_SAMPLE_FILE     } from '../modules/local/createsamplefile/main'
include { BCFTOOLS_REHEADER      } from '../modules/nf-core/bcftools/reheader/main'
include { BCFTOOLS_VIEW          } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_MERGE         } from '../modules/nf-core/bcftools/merge/main'
include { BCFTOOLS_ANNOTATE      } from '../modules/nf-core/bcftools/annotate/main'
include { VCF2COUNTS             } from '../modules/local/vcf2counts/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_vcftocounts_pipeline'

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

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // add index to non-indexed VCFs
    //
    (ch_has_index, ch_has_no_index) = ch_samplesheet
        .map{ it -> [ it[0] + [ name: it[1].simpleName ], it[1], it[2] ] }
        .branch{
            has_index: it[2]
            to_index: !it[2]
                [it[0], it[1]]
        }

    TABIX_TABIX ( ch_has_no_index )

    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    ch_indexed = ch_has_no_index.join(TABIX_TABIX.out.tbi)

    // Join both channels back together
    ch_vcf_tbi = ch_has_index.mix(ch_indexed)

    //
    // Convert gvcfs to vcfs
    //
    (ch_gvcf, ch_normal_vcf) = ch_vcf_tbi.branch{
            gvcf: it[0].gvcf
            vcf: !it[0].gvcf
        }


    GATK4_GENOTYPEGVCFS(
        ch_gvcf.map{ it -> [ it[0], it[1], it[2], [], [] ] },
        fasta.map{ it -> [ [ id:it.baseName ], it ] },
        fai.map{ it -> [ [ id:it.baseName ], it ] },
        dict.map{ it -> [ [ id:it.baseName ], it ] },
        [[],[]], // dbsnp
        [[],[]] // dbsnp_tbi
    )

    ch_vcf_index = GATK4_GENOTYPEGVCFS.out.vcf
            .join(GATK4_GENOTYPEGVCFS.out.tbi)

    ch_vcf = ch_normal_vcf.mix(ch_vcf_index)

    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)

    if (params.filter != null) {
        //
        // Filter VCFs for pattern given in params.filter
        //
        BCFTOOLS_VIEW (
            ch_vcf,
            [], // regions
            [], // targets
            [] // samples
        )

        ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

        ch_filtered_vcf = BCFTOOLS_VIEW.out.vcf
                .join(BCFTOOLS_VIEW.out.tbi)
                .map { meta, vcf, tbi -> [ meta, [ vcf, tbi ] ] }

    } else {
        ch_filtered_vcf = ch_vcf
    }

    //
    // Concatenate converted VCFs if the entries for "id" and "label" are the same
    //
    (ch_single_vcf, ch_multiple_vcf) = ch_filtered_vcf
        .map { meta, vcf, tbi ->
            [ [meta.id, meta.label], meta, vcf, tbi]
        }
        .groupTuple(by: 0)
        .map { _id, metas, vcfs, tbis ->
            def meta = metas[0].findAll { k, _v -> k != 'name' }  // Take the first meta without filename, they should all be the same for a given ID
            [meta, vcfs.flatten(), tbis.flatten()]
        }.branch {
            single:   it[1].size == 1
            multiple: it[1].size > 1
        }

    BCFTOOLS_CONCAT( ch_multiple_vcf )

    ch_vcf_index = BCFTOOLS_CONCAT.out.vcf
            .join(BCFTOOLS_CONCAT.out.tbi)

    ch_vcf_concat = ch_single_vcf
            .mix(ch_vcf_index)

    ch_versions = ch_versions.mix(BCFTOOLS_CONCAT.out.versions)

    if (params.rename) {
        // Create the sample file and add them to ch_vcf_concat
        CREATE_SAMPLE_FILE(ch_vcf_concat.map{ it -> it[0] })

        ch_vcf_sample = ch_vcf_concat
            .join(CREATE_SAMPLE_FILE.out.samplefile)
            .map { meta, vcf, _tbi, samplefile -> [ meta, vcf, [], samplefile ] }

        ch_versions = ch_versions.mix(CREATE_SAMPLE_FILE.out.versions)

        // Rename samples in vcf with the label
        BCFTOOLS_REHEADER(
            ch_vcf_sample,
            [[],[]]
        )

        ch_vcf_index_rh = BCFTOOLS_REHEADER.out.vcf
                .join(BCFTOOLS_REHEADER.out.index)

        ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)
    } else {
        ch_vcf_index_rh = ch_vcf_concat
    }

    //
    // Merge multiple VCFs per sample (patient) with BCFTOOLS_MERGE
    //

    // Bring all vcfs from one sample into a channel
    // Branch based on the number of VCFs per sample
    (ch_single_id, ch_multiple_id) = ch_vcf_index_rh
        .map { meta, vcf, tbi ->
            [meta.id, meta, vcf, tbi]
        }
        .groupTuple(by: 0)
        .map { _id, metas, vcfs, tbis ->
            def meta = metas[0].findAll { k, _v -> k != 'name' }  // Take the first meta without filename, they should all be the same for a given ID
            [meta, vcfs.flatten(), tbis.flatten()]
        }.branch {
            single:   it[1].size == 1
            multiple: it[1].size > 1
        }

    // Run BCFTOOLS_MERGE only on samples with multiple VCFs
    BCFTOOLS_MERGE(
        ch_multiple_id,
        [[],[]], // fasta reference only needed for gvcf
        [[],[]], // fasta.fai reference only needed for gvcf
        [[],[]] // bed
    )

    // Merge the results back into a single channel
    ch_merged_vcfs = ch_single_id.mix(BCFTOOLS_MERGE.out.vcf)

    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    //
    // remove any IDs from the ID column of the VCF
    //
    if (params.removeIDs) {

        BCFTOOLS_ANNOTATE(
            ch_merged_vcfs.map{ it -> [it[0], it[1], [], [], []] },
            [],
            []
        )

        ch_removedIDs_vcfs = ch_single_id.mix(BCFTOOLS_ANNOTATE.out.vcf)

        ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions)
    } else {
        ch_removedIDs_vcfs = ch_merged_vcfs
    }


    //
    // Convert VCFs to Count Matrices
    //
    VCF2COUNTS(
        ch_removedIDs_vcfs.map{ it -> [it[0], it[1]] }
    )

    ch_versions = ch_versions.mix(VCF2COUNTS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'vcftocounts_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    csv            = VCF2COUNTS.out.csv             // channel: *.csv
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
