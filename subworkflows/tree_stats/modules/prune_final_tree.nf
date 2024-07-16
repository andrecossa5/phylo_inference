// PRUNE_FINAL_TREE module
nextflow.enable.dsl = 2 

//

process PRUNE_FINAL_TREE {

    tag "${sample}: ${filtering_key}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering_key}", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        path(tree)
 
    output:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        path("nodes.csv"), 
        path("edges.csv"), emit: final_tree
    
    script:
    """
    python ${baseDir}/bin/process_tree/to_csv.py ${input_folder}
    Rscript ${baseDir}/bin/process_tree/prune_final_tree.r \
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
    touch edges.csv
    touch nodes.csv
    """

}
