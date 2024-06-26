// PROCESS_TREE

// Include here
nextflow.enable.dsl = 2
include { PRUNE_TREE } from "./modules/prune_tree.nf"
include { TREE_STATS } from "./modules/tree_stats.nf"

// 
 
//----------------------------------------------------------------------------//
// PROCESS_TREE subworkflow
//----------------------------------------------------------------------------//

//

workflow PROCESS_TREE {
    
    take:
        ch_input
        ch_tree

    main:

        ch_process = ch_input.combine(ch_tree, by: [0,1])
        PRUNE_TREE(ch_process)
        TREE_STATS(ch_process)

    emit:

        tree = PRUNE_TREE.out.tree
        stats = TREE_STATS.out.stats

}