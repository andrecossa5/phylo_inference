// Re-build and bootstrap tree 

// Include here
nextflow.enable.dsl = 2
include { FINAL_FILTER } from "./modules/final_filter.nf"
include { IQTREE } from "./modules/iqtree.nf"
include { MPBOOT } from "./modules/mpboot.nf"
include { BOOTSTRAP } from "./modules/bootstrap.nf"
include { CASSIOPEIA } from "../raw_tree/modules/cassiopeia.nf"
include { SUPPORT } from "./modules/support.nf"

// 
 
//----------------------------------------------------------------------------//
// FINAL_TREE subworkflow
//----------------------------------------------------------------------------//

// FINAL_TREE subworkflow
workflow FINAL_TREE {
    
    take:
        pruned_tree  

    main: 
        
        // New input_folder equivalent with only cells_assigned to clones with >1% AF for branch-assigned variants. 
        FINAL_FILTER(pruned_tree)
        ch_input = FINAL_FILTER.out.ch_input

        // Final tree building
        if (params.final_solver == "iqtree") {
            tree = IQTREE(ch_input)
        } else if (params.final_solver == "mpboot") {
            tree = MPBOOT(ch_input) 
        } else {
            replicates = Channel.of( 1..params.cassiopeia_n_boot ).concat(Channel.of( "observed" )) 
            BOOTSTRAP(ch_input.combine(replicates))
            CASSIOPEIA(BOOTSTRAP.out.ch_input)
            SUPPORT(CASSIOPEIA.out.tree.groupTuple(by: [0,1]).map{it -> tuple(it[0], it[1], it[4])})
            tree = ch_input.combine(SUPPORT.out.tree, by:[0,1])
        }
 
    emit:

        final_tree = tree

}


 