// phylo_inference aa
nextflow.enable.dsl = 2
include { FILTER_VARIANTS } from "./subworkflows/filter_variants/main"
include { PREPROCESSING } from "./subworkflows/prep_input/main"
include { RAW_TREE } from "./subworkflows/raw_tree/main"
include { FINAL_TREE } from "./subworkflows/final_tree/main"
include { STATS } from "./subworkflows/tree_stats/main"

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
    RAW_TREE(PREPROCESSING.out.input)
    FINAL_TREE(RAW_TREE.out.pruned_tree)
    STATS(FINAL_TREE.out.final_tree)
    STATS.out.stats.view()

}


// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}
