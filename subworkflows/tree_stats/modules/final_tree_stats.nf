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
        path(cell_assignment), 
        path(var_assignment),
        path(final_tree)
 
    output:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        path("final_tree_stats.csv"), emit: stats
    
    script:
    """
    python ${baseDir}/bin/process_tree/final_stats.py \
    --input_folder ${input_folder} \
    --solver ${params.cassiopeia_solver} \
    --trees "${trees}" \
    --n_cores ${task.cpus}
    """

    stub:
    """
    touch final_tree_stats.csv
    """

}
