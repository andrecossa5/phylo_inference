// tree_inference_workflow

// Include here
nextflow.enable.dsl = 2
include { PREP_INPUT } from "./modules/prep_input.nf"
include { BOOTSTRAP } from "./modules/bootstrap.nf"
include { BUILD_TREE } from "./modules/build_tree.nf"

// 
 
//----------------------------------------------------------------------------//
// tree_inference_workflow subworkflow
//----------------------------------------------------------------------------//

workflow tree_inference_workflow {
    
    take:
        ch_samples 

    main:

        PREP_INPUT(ch_samples.combine(params.filtering))
        BOOTSTRAP(
            PREP_INPUT.out.input_folder
            .combine(Channel.of( 1..params.n_boot ))
            .combine(params.boot_strategy)
        )
        // BUILD_TREE(BOOTSTRAP.out.features)

    emit:
        results = BOOTSTRAP.out.bootstrapped_input

}
