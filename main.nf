// phylo_inference pipeline
nextflow.enable.dsl = 2
include { PREPROCESSING } from "./subworkflows/prep_input/main"
include { BUILD_TREE } from "./subworkflows/build_tree/main"
include { EVALUATE_TREE } from "./subworkflows/evaluate_tree/main"

// Samples channel
ch_samples = Channel
    .fromPath("${params.path_data}/*", type:'dir') 
    .map{ it.getName() }

//

//----------------------------------------------------------------------------//
// phylo_inference workflows and entry points
//----------------------------------------------------------------------------//

//

workflow phylo {

    PREPROCESSING(ch_samples)
    BUILD_TREE(PREPROCESSING.out.ch_input)
    EVALUATE_TREE(
        BUILD_TREE.out.ch_tree, 
        PREPROCESSING.out.original_input
    )

}

//

// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}