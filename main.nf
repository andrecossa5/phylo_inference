// phylo_inference, old 
nextflow.enable.dsl = 2
include { FILTER_VARIANTS } from "./subworkflows/filter_variants/main"
include { PREPROCESSING } from "./subworkflows/prep_input/main"
include { CASSIOPEIA } from "./subworkflows/cassiopeia/main"
include { IQTREE } from "./subworkflows/iqtree/main"
include { MPBOOT } from "./subworkflows/mpboot/main"
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


workflow cassiopeia {

    PREPROCESSING(ch_samples)
    CASSIOPEIA(PREPROCESSING.out.input)
    PROCESS_TREE(PREPROCESSING.out.input, CASSIOPEIA.out.tree)
    PROCESS_TREE.out.tree.view()
    PROCESS_TREE.out.stats.view()

}


//


workflow iqtree {

    PREPROCESSING(ch_samples)
    IQTREE(PREPROCESSING.out.input)
    PROCESS_TREE(PREPROCESSING.out.input, IQTREE.out.tree)
    PROCESS_TREE.out.tree.view()
    PROCESS_TREE.out.stats.view()

}


//


workflow mpboot {

    PREPROCESSING(ch_samples)
    MPBOOT(PREPROCESSING.out.input)
    PROCESS_TREE(PREPROCESSING.out.input, MPBOOT.out.tree)
    PROCESS_TREE.out.tree.view()
    PROCESS_TREE.out.stats.view()

}


// Mock
workflow  {
    
    Channel.of(1,2,3,4) | view

}
