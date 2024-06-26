 // PRUNE_RAW_TREE module

nextflow.enable.dsl = 2 

//

process PRUNE_RAW_TREE {

    tag "${sample}: ${filtering_key}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering_key}", mode: 'copy'

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
    Rscript ${baseDir}/bin/process_tree/prune_raw_tree.r # ...
    """

    stub:
    """
    touch cell_assignment.csv
    touch var_assignment.csv
    """

}
