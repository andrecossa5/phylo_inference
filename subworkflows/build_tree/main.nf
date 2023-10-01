// BUILD_TREE

// Include here
nextflow.enable.dsl = 2
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

        if ( ch_input.map{ it[6] in cassiopeia_solvers } ) {
            BUILD_CASSIOPEIA(ch_input)
            tree = BUILD_CASSIOPEIA.out.tree
        }
        else if ( ch_input.map{ it[6] == "raxml"} ) {
            BUILD_RAXML(ch_input)
            tree = BUILD_RAXML.out.tree
        }
        else if ( ch_input.map{ it[6] == "SCITE"} ) {
            BUILD_SCITE(ch_input)
            tree = BUILD_SCITE.out.tree
        }
        else {
            error "Invalid solver: ${solver}"
        }

    emit:
        ch_tree = tree

}
