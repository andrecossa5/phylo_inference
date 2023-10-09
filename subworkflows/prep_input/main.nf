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
        ch_bootstrap = PREP_INPUT.out.input_folder
            .combine(Channel.of( 1..params.n_boot ))
            .combine(params.boot_strategy)
            .map{[ it[0], it[1], it[4], it[3], it[2] ]}
        BOOTSTRAP(ch_bootstrap)
        ch_distance = BOOTSTRAP.out.bootstrapped_input
            .combine(params.distance_based_solver)
            .combine(params.metric)
        ch_others = BOOTSTRAP.out.bootstrapped_input
            .combine(params.other_solver)
            .combine(Channel.of( 'None' ))
            
    emit:
        ch_input = ch_distance
            .concat(ch_others)
            .map{[ it[0], it[1], it[2], it[3], it[5], it[6], it[4] ]}
        original_input = PREP_INPUT.out.input_folder

}
