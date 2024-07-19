// FILTER module

nextflow.enable.dsl = 2

process FILTER_AFM {

    tag "${sample}, filtering_key: ${filtering_key}"

    input:
    tuple val(sample), val(filtering_key)
    
    output:
    tuple val(sample), 
        path("${sample}_${filtering_key}_dataset_df.csv"),
        path("${sample}_${filtering_key}_vars_df.csv"), emit: stats
    tuple val(sample), 
        val(filtering_key), 
        path("${sample}_${filtering_key}_distances.png"),
        path("${sample}_${filtering_key}_tree_raw.png"),
        path("${sample}_${filtering_key}_tree_muts.png"), emit: plots
    
    script:
    """
    python ${baseDir}/bin/filter_variants/get_stats.py \
    --path_data ${params.path_data} \
    --nmads ${params.nmads} \
    --sample_name ${sample} \
    --path_filtering ${params.path_filtering} \
    --filtering_key ${filtering_key} \
    --lineage_column ${params.lineage_column} \
    --solver ${params.cassiopeia_solver} \
    --metric ${params.cassiopeia_metric} \
    --ncores ${task.cpus} \
    --path_priors ${params.path_priors} \
    --path_meta ${params.path_meta} \
    --cell_file ${params.cell_file}
    """

    stub:
    """
    touch "${sample}_${filtering_key}_dataset_df.csv"
    touch "${sample}_${filtering_key}_vars_df.csv"
    touch "${sample}_${filtering_key}_distances.png"
    touch "${sample}_${filtering_key}_tree_raw.png"
    touch "${sample}_${filtering_key}_tree_muts.png"
    """

}
