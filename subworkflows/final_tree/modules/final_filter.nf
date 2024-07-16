// FINAL_FILTER module
nextflow.enable.dsl = 2

process FINAL_FILTER {

    tag "${sample}: ${filtering_key}"

    input:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        path(cell_assignment),
        path(var_assignment)

    output:
    tuple val(sample), 
        val(filtering_key),
        path('filtered_input'), emit: ch_input
    
    script:
    """
    python ${baseDir}/bin/process_tree/final_filter.py \
    --path_input ${input_folder} \
    --filtering_key ${filtering_key} \
    --path_filtering ${params.path_filtering} \
    --path_cell_assignment ${cell_assignment} \
    --path_var_assignment ${var_assignment} \
    --min_assignment_prob ${params.min_assignment_prob}
    """

    stub:
    """
    mkdir filtered_input
    """

}
