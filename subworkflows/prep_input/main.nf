// PREPROCESSING

// Include here
nextflow.enable.dsl = 2
include { PREP_INPUT } from "./modules/prep_input.nf"
include { BOOTSTRAP } from "./modules/bootstrap.nf"

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

        // Input prep: variants filtering
        ch_filtering = Channel.fromPath(params.path_filtering)
                        .map { file -> new groovy.json.JsonSlurper().parse(file).keySet() }
                        .flatMap()
        PREP_INPUT(ch_samples.combine(ch_filtering))

        // Bootstrap filtered character matrices
        ch_bootstrap = PREP_INPUT.out.input_folder
                        .combine(params.boot_strategy)
                        .combine(Channel.of( 1..params.n_boot ))
        BOOTSTRAP(ch_bootstrap)
            
    emit:
        boot_input = BOOTSTRAP.out.bootstrapped_input.combine(params.solver)
        original_input = PREP_INPUT.out.input_folder

}
