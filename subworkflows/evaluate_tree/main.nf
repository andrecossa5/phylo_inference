// EVALUATE_TREE

// Include here
nextflow.enable.dsl = 2
include { INTERNAL_CONSISTENCY } from "./modules/internal_consistency.nf"
include { SOLVERS_CONSISTENCY } from "./modules/solvers_consistency.nf"

// 
 
//----------------------------------------------------------------------------//
// EVALUATE_TREE subworkflow
//----------------------------------------------------------------------------//

//

workflow EVALUATE_TREE {
    
    take:
        ch_tree

    main:

        // Internal consistency within boot trees
        INTERNAL_CONSISTENCY(ch_tree.groupTuple(by: [0,1,3,5]))
        // Concordance (same sample and input features, different solvers)
        SOLVERS_CONSISTENCY(
            ch_tree
            .filter { it[6].endsWith('rep_observed.newick') }
            .groupTuple(by: [0,1])
        )

    emit:
    
        results = SOLVERS_CONSISTENCY.out.results


}