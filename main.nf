// phylo_inference, old 
nextflow.enable.dsl = 2
include { FILTER_VARIANTS } from "./subworkflows/filter_variants/main"
include { PREPROCESSING } from "./subworkflows/prep_input/main"
include { CASSIOPEIA } from "./subworkflows/cassiopeia/main"
// include { IQTREE } from "./subworkflows/iqtree/main"
// include { MPBOOT } from "./subworkflows/mpboot/main"
include { PROCESS_TREE } from "./subworkflows/process_tree/main"

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
    CASSIOPEIA(PREPROCESSING.out.input)
    PROCESS_TREE(PREPROCESSING.out.input, CASSIOPEIA.out.tree)

    if (params.) {
        println 'Go on, condition matched.'
    } else {
        println 'Condition NOT matched. Skip execution.'
    }

}


// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}
