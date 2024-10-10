// PREP_MAESTER module

nextflow.enable.dsl = 2 

//

process ONESAMPLE {

    tag "${sample}: ${job_id}"
    publish "${params.outdir}/${sample}"

    input:
    tuple val(job_id), 
        val(sample), 
        path(path_afm), 
        val(char_filtering_key),
        val(cell_filtering_key),
        val(bin_key),
        val(tree_key),
        val(cell_file)
 
    output:
    tuple val(job_id), 
        val(sample), 
        path("${job_id}_stats.pickle"), emit: stats
    
    script:
    """
    python ${baseDir}/bin/pp/onesample.py \
    --path_afm ${path_afm} \
    --job_id ${job_id} \
    --cell_filter ${params.cell_filter} \
    --min_cell_number ${params.min_cell_number} \
    --min_cov ${params.min_cov} \
    --min_var_quality ${params.min_var_quality} \
    --min_frac_negative ${params.min_frac_negative} \
    --min_n_positive ${min_n_positive} \
    --af_confident_detection ${af_confident_detection} \
    --min_n_confidently_detected ${min_n_confidently_detected} \
    --min_mean_AD_in_positives ${min_mean_AD_in_positives} \
    --min_mean_DP_in_positives ${params.min_mean_DP_in_positives} \
    --t_prob ${params.t_prob} \
    --t_vanilla ${params.t_vanilla} \
    --min_AD ${min_AD} \
    --min_cell_prevalence ${params.min_cell_prevalence} \
    --bin_method ${bin_method} \
    --solver ${params.solver} \
    --lineage_column ${params.lineage_column} \
    --ncores ${task.cpus} \
    --cell_file ${params.cell_file} \
    --path_dbSNP ${params.path_dbSNP} \
    --path_REDIdb ${params.}
    """

    stub:
    """
    ${job_id}_stats.pickle
    """

}
