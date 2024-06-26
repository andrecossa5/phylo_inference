// FILTER_VARIANTS

// Include here
nextflow.enable.dsl = 2
include { FILTER_AFM } from "./modules/filter_afm.nf"
include { publish_images } from "./modules/publish_plots.nf"
include { collapse_output } from "./modules/collapse_output.nf"

// 
 
//----------------------------------------------------------------------------//
// FILTER_VARIANTS subworkflow
//----------------------------------------------------------------------------//

workflow FILTER_VARIANTS {
    
    take:
        ch_samples  

    main:

        // Variants filtering
        ch_filtering = Channel.fromPath(params.path_filtering)
                        .map { file -> new groovy.json.JsonSlurper().parse(file).keySet() }
                        .flatMap()
        FILTER_AFM(ch_samples.combine(ch_filtering))

        // Collapse outputs
        ch_grouped = FILTER_AFM.out.stats
            .groupTuple(by: 0)
            .flatMap { sample, datasetPaths, varsPaths ->
                [
                    [sample, 'dataset', datasetPaths],
                    [sample, 'vars', varsPaths]
                ]
            }
        collapse_output(ch_grouped)

        // Publish plots
        publish_images(FILTER_AFM.out.plots)
        
            
    emit:
        results = collapse_output.out.csv

}