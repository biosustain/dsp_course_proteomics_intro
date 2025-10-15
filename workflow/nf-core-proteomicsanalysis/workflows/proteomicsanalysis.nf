/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_proteomicsanalysis_pipeline'

include { VUEGEN } from '../modules/nf-core/vuegen/main'   
include { ACORE } from '../modules/local/acore/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROTEOMICSANALYSIS {

    // take:
    // ch_input // pass on output channels from other process/subworkflows/modules

    main:

    ch_versions = Channel.empty()
    
    // MODULE: Run Analysis notebook to generate result folder
    
    input_csv_ch = Channel.value(params.input_csv)
    nb_ch = Channel.value(params.nb)

    ACORE(
        input_csv_ch,
        nb_ch
    )

    report_type_ch = Channel.value(params.report_type)
    input_type_ch = Channel.value(params.input_type) // 'directory', 'config', etc.

    VUEGEN(
        input_type_ch,    
        ACORE.out.report_files, // Use the results folder created by Acore module
        report_type_ch,
    )
    //
    // Collate and save software versions
    //
    ch_versions = ch_versions.mix(VUEGEN.out.versions.first())
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'proteomicsanalysis_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )


    emit: versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
