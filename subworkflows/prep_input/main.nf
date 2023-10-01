// PREPROCESSING

// Include here
nextflow.enable.dsl = 2
include { PREP_INPUT } from "./modules/prep_input.nf"
include { BOOTSTRAP } from "./modules/bootstrap.nf"

// 
 
//----------------------------------------------------------------------------//
// PREPROCESSING subworkflow
//----------------------------------------------------------------------------//

// PREPROCESSING subworkflow
workflow PREPROCESSING {
    
    take:
        ch_samples  

    main:

        PREP_INPUT(ch_samples.combine(params.filtering))
        BOOTSTRAP(
            PREP_INPUT.out.input_folder
            .combine(Channel.of( 1..params.n_boot ))
            .combine(params.boot_strategy)
        )
        ch_distance = BOOTSTRAP.out.bootstrapped_input
            .combine(params.distance_based_solver)
            .combine(params.metric)
        ch_others = BOOTSTRAP.out.bootstrapped_input
            .combine(params.other_solver)
            .combine(Channel.of( 'None' ))
            
    emit:
        ch_input = ch_distance.concat(ch_others)
        original_input = PREP_INPUT.out.input_folder

}
