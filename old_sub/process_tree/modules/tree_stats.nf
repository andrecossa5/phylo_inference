// TREE_STATS module

nextflow.enable.dsl = 2 

//

process TREE_STATS {

    tag "${sample}: ${filtering_key}, ${boot_method}, ${solver}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering_key}/${solver}/", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering_key),
        path(input_folder),
        val(solver),
        path(tree)

    output:
    tuple val(sample), 
        val(filtering_key), 
        val(solver),
        path("${filtering_key}_${solver}_stats.csv"), emit: stats
    
    script:
    """
    python ${baseDir}/bin/process_tree/tree_stats.py \
    --path_input ${input_folder} \
    --path_meta ${params.path_meta} \
    --path_priors ${params.path_priors} \
    --path_filtering ${params.path_filtering} \
    --filtering_key ${filtering_key} \
    --solver ${solver} \
    --sample ${sample} \
    --n_cores ${task.cpus} \
    --lineage_column ${params.lineage_column}
    """

    stub:
    """ 
    touch ${filtering_key}_${solver}_stats.csv
    """

}
