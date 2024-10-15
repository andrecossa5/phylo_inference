// ONESAMPLE module

nextflow.enable.dsl = 2 

//

process ONESAMPLE {

    tag "${sample}: tuning${job_id}"
    publishDir "${params.outdir}/${sample}", mode: "copy"

    input:
    tuple val(sample), 
        path(path_afm), 
        val(min_n_positive),
        val(af_confident_detection),
        val(min_n_confidently_detected),
        val(min_mean_AD_in_positives),
        val(min_AD),
        val(bin_method),
        val(t_prob),
        val(min_cell_prevalence),
        val(job_id)

 
    output:
    tuple val(job_id), 
        val(sample), 
        path("tuning${job_id}_stats.pickle"), emit: stats
    
    script:
    """
    python ${baseDir}/bin/pp/onesample.py \
    --path_afm ${path_afm} \
    --job_id ${job_id} \
    --cell_filter ${params.cell_filter} \
    --filtering ${params.filtering} \
    --min_cell_number ${params.min_cell_number} \
    --min_cov ${params.min_cov} \
    --min_var_quality ${params.min_var_quality} \
    --min_frac_negative ${params.min_frac_negative} \
    --min_n_positive ${min_n_positive} \
    --af_confident_detection ${af_confident_detection} \
    --min_n_confidently_detected ${min_n_confidently_detected} \
    --min_mean_AD_in_positives ${min_mean_AD_in_positives} \
    --min_mean_DP_in_positives ${params.min_mean_DP_in_positives} \
    --t_prob ${t_prob} \
    --t_vanilla ${params.t_vanilla} \
    --min_AD ${min_AD} \
    --min_cell_prevalence ${min_cell_prevalence} \
    --bin_method ${bin_method} \
    --solver ${params.solver} \
    --lineage_column ${params.lineage_column} \
    --ncores ${task.cpus} \
    --cell_file ${params.cell_file} \
    --path_dbSNP ${params.path_dbSNP} \
    --path_REDIdb ${params.path_REDIdb}
    """

    stub:
    """
    touch tuning${job_id}_stats.pickle
    """

}