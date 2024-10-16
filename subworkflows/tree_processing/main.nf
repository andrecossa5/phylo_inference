// Build a cell phylogeny

// Include here
nextflow.enable.dsl = 2
include { PROCESS_TREE } from "./modules/process_tree.nf"
include { TREE_METRICS } from "./modules/tree_metrics.nf"

// 
 
//----------------------------------------------------------------------------//
// tree subworkflow
//----------------------------------------------------------------------------//

//

workflow process_tree {
    
    take: 
        ch_input 
        ch_tree  

    main: 

        ch_annot = ch_input.flatMap { 
                job_id, sample, key1, key2, input_folder, replicates, distances ->
                replicates.indices.collect { i ->
                    tuple(job_id, sample, input_folder, replicates[i], distances[i])
                }
            }
            .filter(it -> it[3]=="observed")
            .combine(ch_tree, by:[0,1])

        PROCESS_TREE(ch_annot)
        TREE_METRICS(PROCESS_TREE.out.annotated_tree)

    emit:

        metrics = TREE_METRICS.out.metrics

}


 