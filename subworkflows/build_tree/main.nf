// BUILD_TREE

// Include here
nextflow.enable.dsl = 2
include { CALC_DISTANCES } from "./modules/calculate_distances.nf"
include { BUILD_CASSIOPEIA } from "./modules/build_cassiopeia.nf"
include { BUILD_RAXML } from "./modules/build_raxml.nf"
include { BUILD_SCITE } from "./modules/build_scite.nf"

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

        cassiopeia_solvers = [
            "NJ", "UPMGA", "max_cut", "greedy", 
            "spectral", "shared_muts"
        ]

        // Distances
        ch_dist = CALC_DISTANCES(
            ch_input
            .map{[ it[0], it[1], it[2], it[3], it[5], it[6] ]}
            .unique()
        )

        // Bootstrapped trees
        if ( ch_input.map{ it[4] in cassiopeia_solvers } ) {
            BUILD_CASSIOPEIA(
                ch_input
                .map{[ it[0], it[1], it[2], it[3], it[5], it[4], it[6] ]}
                .combine(CALC_DISTANCES.out.dist, by:[0,1,2,3,4])
            )
            tree = BUILD_CASSIOPEIA.out.tree
        }
        else if ( ch_input.map{ it[4] == "raxml"} ) {
            BUILD_RAXML(ch_input)
            tree = BUILD_RAXML.out.tree
        }
        else if ( ch_input.map{ it[4] == "SCITE"} ) {
            BUILD_SCITE(ch_input)
            tree = BUILD_SCITE.out.tree
        }
        else {
            error "Invalid solver: ${solver}"
        }

    emit:
        ch_tree = tree

}