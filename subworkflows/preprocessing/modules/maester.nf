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
        val(bin_key), 
        val(tree_key), 
        path("${job_id}_pp"), emit: input_folder
    
    script:
    """
    python ${baseDir}/bin/pp/MAESTER.py \
    --path_afm ${ch_matrix} \
    --path_meta ${params.path_meta} \
    --path_char_filtering ${params.path_char_filtering} \
    --path_cell_filtering ${params.path_cell_filtering} \
    --path_bin ${params.path_bin} \
    --path_tree ${params.path_distance_tree} \
    --char_filtering_key ${char_filtering_key} \
    --cell_filtering_key ${cell_filtering_key} \
    --bin_key ${bin_key} \
    --tree_key ${tree_key} \
    --sample ${sample} \
    --job_id ${job_id} \
    --lineage_column ${params.lineage_column} \
    --n_cores ${task.cpus} \
    --cell_file ${cell_file} \
    --path_dbSNP ${params.path_dbSNP} \
    --path_REDIdb ${params.path_REDIdb}
    """

    stub:
    """
    mkdir ${job_id}_pp
    cd ${job_id}_pp
    touch AD.npz
    touch DP.npz
    touch bin_ops.csv
    touch cell_filtering_ops.csv
    touch cell_meta.csv
    touch char_filtering_ops.csv
    touch char_meta.csv
    touch dataset_df.csv
    touch tree_ops.csv
    touch genotypes.fa
    """

}
