// phylo_inference 
nextflow.enable.dsl = 2
include { FILTER_VARIANTS } from "./subworkflows/filter_variants/main"
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


workflow muts {

    FILTER_VARIANTS(ch_samples)
    FILTER_VARIANTS.out.results.view()

}


//


workflow phylo {

    PREPROCESSING(ch_samples)
    BUILD_TREE(PREPROCESSING.out.input)
    EVALUATE_TREE(BUILD_TREE.out.ch_tree)
    // EVALUATE_TREE.out.results.view()

}


// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}
