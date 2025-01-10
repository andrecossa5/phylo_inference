// MAESTER and MAESTER_TUNE module

nextflow.enable.dsl = 2 

//

process MAESTER_TUNE {

    tag "${sample}: ${job_id}"

    input:
    tuple val(job_id), val(sample), val(ch_matrix)
 
    output:
    tuple val(job_id), val(sample), path("afm.h5ad"), emit: afm
    
    script:
    """
    python ${baseDir}/bin/pp/MAESTER.py \
    --path_afm ${ch_matrix} \
    --path_pickles ${params.path_pickles} \
    --sample ${sample} \
    --job_id ${job_id} \
    --n_cores ${task.cpus} \
    --path_dbSNP ${params.path_dbSNP} \
    --path_REDIdb ${params.path_REDIdb}
    """

    stub:
    """
    touch afm.h5ad
    """

}


//


process MAESTER {

    tag "${sample}: ${job_id}"

    input:
    tuple val(job_id), val(sample), val(ch_matrix)
 
    output:
    tuple val(job_id), val(sample), path("afm.h5ad"), emit: afm
    
    script:
    """
    python ${baseDir}/bin/pp/MAESTER.py \
    --path_afm ${ch_matrix} \
    --sample ${sample} \
    --cell_filter ${params.cell_filter} \
    --filtering ${params.filtering} \
    --min_cell_number ${params.min_cell_number} \
    --min_cov ${params.min_cov} \
    --min_var_quality ${params.min_var_quality} \
    --min_frac_negative ${params.min_frac_negative} \
    --min_n_positive ${params.min_n_positive} \
    --af_confident_detection ${params.af_confident_detection} \
    --min_n_confidently_detected ${params.min_n_confidently_detected} \
    --min_mean_AD_in_positives ${params.min_mean_AD_in_positives} \
    --min_mean_DP_in_positives ${params.min_mean_DP_in_positives} \
    --t_prob ${params.t_prob} \
    --t_vanilla ${params.t_vanilla} \
    --min_AD ${params.min_AD} \
    --min_cell_prevalence ${params.min_cell_prevalence} \
    --bin_method ${params.bin_method} \
    --lineage_column ${params.lineage_column} \
    --solver ${params.cassiopeia_solver} \
    --metric ${params.distance_metric} \
    --n_cores ${task.cpus} \
    --path_dbSNP ${params.path_dbSNP} \
    --path_REDIdb ${params.path_REDIdb}
    """

    stub:
    """
    touch afm.h5ad
    """

}
