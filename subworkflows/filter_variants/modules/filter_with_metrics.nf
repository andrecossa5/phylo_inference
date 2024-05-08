// FILTER module

nextflow.enable.dsl = 2

process FILTER {

    tag "${sample}:job_${id}"

    input:
    tuple val(id),
        val(sample), 
        val(min_site_cov),
        val(min_var_quality),
        val(min_frac_negative),
        val(min_n_positive),
        val(low_confidence_af),
        val(high_confidence_af),
        val(min_prevalence_low_confidence_af),
        val(min_cells_high_confidence_af)
    
    output:
    tuple val(sample), 
        path('dataset_df.csv'),
        path('vars_df.csv'), emit: stats
    
    script:
    """
    python ${baseDir}/bin/filter_variants/get_stats.py \
    --path_data ${params.path_data} \
    --min_site_cov ${min_site_cov} \
    --min_var_quality ${min_var_quality} \
    --min_frac_negative ${min_frac_negative} \
    --min_n_positive ${min_n_positive} \
    --low_confidence_af ${low_confidence_af} \
    --high_confidence_af ${high_confidence_af} \
    --min_prevalence_low_confidence_af ${min_prevalence_low_confidence_af} \
    --min_cells_high_confidence_af ${min_cells_high_confidence_af} \
    --lineage_column ${params.lineage_column} \
    --solver ${params.solver} \
    --metric ${params.metric} \
    --ncores ${task.ncpus} \
    --path_priors ${params.path_priors} \
    --path_meta ${params.path_meta} \
    --job_id ${id}
    """

    stub:
    """
    touch job_df.csv
    touch dataset_df.csv
    touch vars_df.csv
    """

}
