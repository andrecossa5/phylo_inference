// Process the final tree

// Include here
nextflow.enable.dsl = 2
include { PRUNE_FINAL_TREE } from "./modules/prune_final_tree.nf"
include { FINAL_TREE_STATS } from "./modules/final_tree_stats.nf"
include { TEST_LINEAGE_ASSOCIATION } from "./modules/test_lineage_association.nf"

// 
 
//----------------------------------------------------------------------------//
// STATS subworkflow
//----------------------------------------------------------------------------//

//

// STATS subworkflow
workflow STATS {
    
    take: 
        final_tree   

    main: 
        
        PRUNE_FINAL_TREE(final_tree)
        FINAL_TREE_STATS(PRUNE_FINAL_TREE.out.final_tree)

        // Test association with GT lineage covariate, if necessary
        if (params.test_association_with_lineage) {
            TEST_LINEAGE_ASSOCIATION(PRUNE_FINAL_TREE.out.final_tree)
        }

    emit:

        stats = FINAL_TREE_STATS.out.stats

}


 

