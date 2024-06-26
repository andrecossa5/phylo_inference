// Build the first tree draft

// Include here
nextflow.enable.dsl = 2
include { CASSIOPEIA } from "./modules/cassiopeia.nf"
include { PRUNE_RAW_TREE } from "./modules/prune_raw_tree.nf"

// 
 
//----------------------------------------------------------------------------//
// CASSIOPEIA subworkflow
//----------------------------------------------------------------------------//

//

// RAW_TREE subworkflow
workflow RAW_TREE {
    
    take: 
        ch_input   

    main: 
        
        CASSIOPEIA(ch_input.combine(Channel.of( "raw" )))
        PRUNE_RAW_TREE(CASSIOPEIA.out.tree)

    emit:

        pruned_tree = PRUNE_RAW_TREE.out.pruned_tree

}


 