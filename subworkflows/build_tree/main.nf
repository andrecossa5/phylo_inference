// BUILD_TREE

// Include here
nextflow.enable.dsl = 2
include { CALC_DISTANCES } from "./modules/calculate_distances.nf"
include { BUILD_CASSIOPEIA } from "./modules/build_cassiopeia.nf"
include { BUILD_MPBOOT } from "./modules/build_mpboot.nf"
include { BUILD_IQTREE } from "./modules/build_iqtree.nf"

// 
 
//----------------------------------------------------------------------------//
// BUILD_TREE subworkflow
//----------------------------------------------------------------------------//

//

// BUILD_TREE subworkflow
workflow BUILD_TREE {
    
    take:
        ch_input  

    main:

        // Cassiopeia solvers
        cassiopeia_solvers = ["NJ", "UPMGA", "max_cut", "greedy", "spectral", "shared_muts"]

        // Choose tree building method and construct trees
        if ( ch_input.map{ it[5] in cassiopeia_solvers } ) {

            CALC_DISTANCES(ch_input)
            BUILD_CASSIOPEIA(CALC_DISTANCES.out.dist)
            ch_tree = BUILD_CASSIOPEIA.out.tree

        }
        else if ( ch_input.map{ it[5] == "iqtree"} ) {
            BUILD_IQTREE(ch_input)
            ch_tree = BUILD_IQTREE.out.tree
        }
        else if ( ch_input.map{ it[5] == "MPBoot"} ) {
            BUILD_MPBOOT(ch_input)
            ch_tree = BUILD_MPBOOT.out.tree
        }
        else {
            error "Invalid solver: ${solver}"
        }

    emit:

        ch_tree = ch_tree

}