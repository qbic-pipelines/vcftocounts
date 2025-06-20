// Only run when --filter or --subset is set
workflow WRITE_META_TO_CSV {
    take:
        ch_vcf                  // channel: [mandatory] meta, vcf, tbi
        csv_subfolder           //
        outdir                  //

    main:
        // Creating csv files to restart from this step
        ch_vcf .collectFile(keepHeader: true, skip: 1, sort: true, storeDir: "${outdir}/csv") { meta, file, index ->
            def sample         = meta.id
            def label          = meta.label
            def numVarBefore   = meta.numVarBefore
            def numVarAfter    = meta.numVarAfter
            def suffix_aligned = save_output_as_bam ? "bam" : "cram"
            def suffix_index   = save_output_as_bam ? "bam.bai" : "cram.crai"
            file   = "${outdir}/preprocessing/${csv_subfolder}/${sample}/${file.baseName}.vcf.gz"
            index   = "${outdir}/preprocessing/${csv_subfolder}/${sample}/${index.baseName.minus(".vcf.gz")}.tbi"

            ["markduplicates_no_table.csv", "patient,sex,status,sample,${type},${type_index}\n${patient},${sex},${status},${sample},${file},${index}\n"]
        }
}
