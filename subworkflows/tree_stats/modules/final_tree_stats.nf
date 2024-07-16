// FINAL_TREE_STATS module
nextflow.enable.dsl = 2 

//

process FINAL_TREE_STATS {

    tag "${sample}: ${filtering_key}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering_key}", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        path(nodes), 
        path(edges)
 
    output:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        path("final_tree_stats.csv"),
        path("annotated_tree.pickle"), emit: stats
    
    script:
    """
    python ${baseDir}/bin/process_tree/final_stats.py \
    --path_input ${input_folder} \
    --path_nodes ${nodes} \
    --path_edges ${edges} \
    --metric ${params.distance_metric} \
    --ncores ${task.cpus}
    """

    stub:
    """
    touch final_tree_stats.csv
    touch annotated_tree.pickle
    """

}
