 // PRUNE_RAW_TREE module

nextflow.enable.dsl = 2 

//

process PRUNE_RAW_TREE {

    tag "${sample}: ${filtering_key}"

    input:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        val(tree_name),
        path(tree)
 
    output:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        path("cell_assignment.csv"), 
        path("var_assignment.csv"), emit: pruned_tree
    
    script:
    """
    python ${baseDir}/bin/process_tree/to_csv.py ${input_folder}
    Rscript ${baseDir}/bin/process_tree/prune_raw_tree.r \
    --tree ${tree} \
    --AD AD.csv \
    --DP DP.csv \
    --af_t 0.05 \
    --ncores ${task.cpus} \
    --prob_cut 0.3 \
    --min_cell_clone ${params.min_cell_clone}
    """

    stub:
    """
    touch cell_assignment.csv
    touch var_assignment.csv
    """

}
