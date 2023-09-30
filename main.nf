// phylo_inference pipeline
nextflow.enable.dsl = 2
include { PREPROCESSING } from "./subworkflows/prep_input/main"
include { BUILD_TREE } from "./subworkflows/build_tree/main"

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
    BUILD_TREE.out.reconstructed_tree.view()

}

//

// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}