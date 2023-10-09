// EVALUATE_TREE

// Include here
nextflow.enable.dsl = 2
include { CALC_DISTANCES } from "./modules/calculate_distances.nf"
include { BUILD_OBSERVED_CASSIOPEIA } from "./modules/build_obs_cassiopeia.nf"
include { BUILD_OBSERVED_RAXML } from "./modules/build_obs_raxml.nf"
include { BUILD_OBSERVED_SCITE } from "./modules/build_obs_scite.nf"
include { EVALUATE_I } from "./modules/evaluate_I.nf"
include { EVALUATE_II } from "./modules/evaluate_II.nf"
include { EVALUATE_III } from "./modules/evaluate_III.nf"

// 
 
//----------------------------------------------------------------------------//
// EVALUATE_TREE subworkflow
//----------------------------------------------------------------------------//

//

workflow EVALUATE_TREE {
    
    take:
        ch_tree  
        original_input

    main:

        // Channels manipulation
        ch_trees_by_run =  ch_tree.groupTuple(by: [0,1,2,4,5])
        ch_obs_trees = ch_trees_by_run
            .map{[ it[0], it[1], it[4], it[5] ]}
            .unique()
            .combine(original_input, by:[0,1])

        // Build observed trees
        cassiopeia_solvers = [
            "NJ", "UPMGA", "max_cut", "greedy", 
            "spectral", "shared_muts"
        ]

        // Distances
        CALC_DISTANCES(ch_obs_trees.map{[ it[0], it[1], it[3], it[4] ]}.unique())

        // Observed trees
        if ( ch_obs_trees.map{ it[2] in cassiopeia_solvers } ) {
            BUILD_OBSERVED_CASSIOPEIA(
                ch_obs_trees
                .map{[ it[0], it[1], it[3], it[2] ]}
                .combine(CALC_DISTANCES.out.dist, by:[0,1,2])
                .map{[ it[0], it[1], it[3], it[2], it[4], it[5] ]}
            )
            obs_tree = BUILD_OBSERVED_CASSIOPEIA.out.tree
        }
        else if ( ch_obs_trees.map{ it[3] == "raxml"} ) {
            BUILD_RAXML(ch_obs_trees)
            obs_tree = BUILD_OBSERVED_RAXML.out.tree
        }
        else if ( ch_obs_trees.map{ it[3] == "SCITE"} ) {
            BUILD_OBSERVED_SCITE(ch_obs_trees)
            obs_tree = BUILD_OBSERVED_SCITE.out.tree
        }
        else {
            error "Invalid solver: ${solver}"
        }

        // Evaluation

        // Internal consistency
        EVALUATE_I(
            ch_trees_by_run
            .map{[ it[0], it[1], it[4], it[5], it[2], it[6] ]}
            .combine(obs_tree, by: [0,1,2,3])
        )
        // Concordance (same sample and input features)
        EVALUATE_II(obs_tree.groupTuple(by: [0,1]))
        // Evaluate external associations
        EVALUATE_III(obs_tree.combine(original_input, by: [0,1]))

    emit:
        results = obs_tree.combine(original_input, by: [0,1])


}