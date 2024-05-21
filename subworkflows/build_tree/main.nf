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
        
        // Run tree building
        ch_input_tree = ch_input
            .branch {
                cassiopeia: it[5] in cassiopeia_solvers
                iqtree: it[5] == "iqtree"
                mpboot: it[5] == "MPBoot"
            }  
        CALC_DISTANCES(ch_input_tree.cassiopeia)
        BUILD_CASSIOPEIA(CALC_DISTANCES.out.dist)
        BUILD_IQTREE(ch_input_tree.iqtree)
        BUILD_MPBOOT(ch_input_tree.mpboot)

        // Concat channels
        ch_tree = BUILD_CASSIOPEIA.out.tree.concat(BUILD_IQTREE.out.tree, BUILD_MPBOOT.out.tree)

    emit:

        ch_tree = ch_tree

}


