// PREP_MAESTER module

nextflow.enable.dsl = 2 

//

process MAESTER {

    tag "${sample}: ${job_id}"

    input:
    tuple val(job_id), 
        val(sample), 
        val(ch_matrix), 
        val(char_filtering_key),
        val(cell_filtering_key),
        val(bin_key),
        val(tree_key),
        val(cell_file)
 
    output:
    tuple val(job_id), 
        val(sample), 
        path("afm.h5ad"), emit: afm
    
    script:
    """
    python ${baseDir}/bin/pp/MAESTER.py \
    --path_afm ${ch_matrix} \
    --path_char_filtering ${params.path_char_filtering} \
    --path_cell_filtering ${params.path_cell_filtering} \
    --path_bin ${params.path_bin} \
    --path_tree ${params.path_distance_tree} \
    --char_filtering_key ${char_filtering_key} \
    --cell_filtering_key ${cell_filtering_key} \
    --bin_key ${bin_key} \
    --tree_key ${tree_key} \
    --job_id ${job_id} \
    --lineage_column ${params.lineage_column} \
    --n_cores ${task.cpus} \
    --cell_file ${cell_file} \
    --path_dbSNP ${params.path_dbSNP} \
    --path_REDIdb ${params.path_REDIdb}
    """

    stub:
    """
    touch afm.h5ad
    """

}
