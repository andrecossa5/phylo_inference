// BUILD_TREE

// Include here
nextflow.enable.dsl = 2
include { CALC_DISTANCES } from "./modules/calculate_distances.nf"
include { BUILD_CASSIOPEIA } from "./modules/build_cassiopeia.nf"
include { BUILD_RAXML } from "./modules/build_raxml.nf"

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
        // else if ( ch_input.map{ it[5] == "raxml"} ) {
        //     BUILD_RAXML(ch_input)~
        //     tree = BUILD_RAXML.out.tree
        // }
        // else if ( ch_input.map{ it[5] == "SCITE"} ) {
        //     BUILD_SCITE(ch_input)
        //     tree = BUILD_SCITE.out.tree
        // }
        // else {
        //     error "Invalid solver: ${solver}"
        // }

    emit:

        ch_tree = ch_tree

}