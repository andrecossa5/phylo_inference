// build_tree

// Include here
nextflow.enable.dsl = 2
include { CASSIOPEIA } from "./modules/cassiopeia.nf"
include { IQTREE } from "./modules/iqtree.nf"
include { MPBOOT } from "./modules/mpboot.nf"
include { BOOSTER } from "./modules/booster.nf"

// 
 
//----------------------------------------------------------------------------//
// build_tree subworkflow
//----------------------------------------------------------------------------//

//

workflow build_tree {
    
    take: 
        ch_input   

    main: 

        ch_flattened = ch_input.flatMap { 
            job_id, sample, replicates, afms ->
                replicates.indices.collect { i ->
                tuple(job_id, sample, replicates[i], afms[i])
           }
        } 

        if (params.tree_algorithm == "cassiopeia") {
            CASSIOPEIA(ch_flattened)
            trees = CASSIOPEIA.out.tree.groupTuple(by: [0,1])
        } else if (params.tree_algorithm == "mpboot") {
            trees = MPBOOT(ch_flattened.filter(it -> it[2]=="observed"))
        } else if (params.tree_algorithm == "iqtree") {
            trees = IQTREE(ch_flattened.filter(it -> it[2]=="observed"))
        } else {
            println('Provide valid tracing system option! (e.g., cassiopeia, mpboot, iqtree)')
        }

        BOOSTER(trees)

    emit:

        final_tree = BOOSTER.out.final_tree
 
}
 

 