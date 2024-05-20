// FILTER module

nextflow.enable.dsl = 2

process FILTER_AFM {

    tag "${sample}, filtering_key: ${filtering_key}"

    input:
    tuple val(sample), val(filtering_key)
    
    output:
    tuple val(sample), 
        path("${sample}_${filtering_key}_job_df.csv"),
        path("${sample}_${filtering_key}_dataset_df.csv"),
        path("${sample}_${filtering_key}_vars_df.csv"), emit: stats
    
    script:
    """
    python ${baseDir}/bin/filter_variants/get_stats.py \
    --path_data ${params.path_data} \
    --sample_name ${sample} \
    --path_filtering ${params.path_filtering} \
    --filtering_key ${filtering_key}
    --lineage_column ${params.lineage_column} \
    --solver ${params.solver} \
    --metric ${params.metric} \
    --spatial_metrics ${params.spatial_metrics} \
    --ncores ${task.cpus} \
    --path_priors ${params.path_priors} \
    --path_meta ${params.path_meta}
    """

    stub:
    """
    touch "${sample}_${filtering_key}_job_df.csv"
    touch "${sample}_${filtering_key}_dataset_df.csv"
    touch "${sample}_${filtering_key}_vars_df.csv"
    """

}
