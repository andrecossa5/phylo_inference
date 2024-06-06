// PREPROCESSING

// Include here
nextflow.enable.dsl = 2
include { PREP_INPUT } from "./modules/prep_input.nf"

// Import
import groovy.json.JsonSlurper

// 
 
//----------------------------------------------------------------------------//
// PREPROCESSING subworkflow
//----------------------------------------------------------------------------//

// PREPROCESSING subworkflow
workflow PREPROCESSING {
    
    take:
        ch_samples  

    main:

        ch_filtering = Channel.fromPath(params.path_filtering)
                        .map { file -> new groovy.json.JsonSlurper().parse(file).keySet() }
                        .flatMap()
        PREP_INPUT(ch_samples.combine(ch_filtering))
            
    emit:

        input = PREP_INPUT.out.input_folder

} 
