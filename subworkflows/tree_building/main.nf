// build_tree

// Include here
nextflow.enable.dsl = 2
include { CASSIOPEIA } from "./modules/cassiopeia.nf"
include { SUPPORT } from "./modules/support.nf"
include { IQTREE } from "./modules/iqtree.nf"
include { MPBOOT } from "./modules/mpboot.nf"

// 
 
//----------------------------------------------------------------------------//
// build_tree subworkflow
//----------------------------------------------------------------------------//

//

workflow build_tree {
    
    take: 
        ch_input   

    main: 

        if (params.tree_algorithm == "cassiopeia") {
            ch_flattened = ch_input.flatMap { 
                job_id, sample, key1, key2, input_folder, replicates, distances ->
                replicates.indices.collect { i ->
                    tuple(job_id, sample, key1, key2, input_folder, replicates[i], distances[i])
                }
            }
            CASSIOPEIA(ch_flattened)
            final_tree = SUPPORT(CASSIOPEIA.out.tree.groupTuple(by: [0,1,2,3,4]))
        } else if (params.tree_algorithm == "mpboot") {
            final_tree = MPBOOT(ch_input.map{ it -> it[0], it[1], it[4] })
        } else if (params.tree_algorithm == "iqtree") {
            final_tree = IQTREE(ch_input.map{ it -> it[0], it[1], it[4] })
        } else {
            println('Provide valid tracing system option! (e.g., cassiopeia, mpboot, iqtree)')
        }

    emit:

        final_tree = final_tree

}


 