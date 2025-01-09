// PREP_MAESTER module

nextflow.enable.dsl = 2 

//

process MAESTER {

    tag "${sample}: ${job_id}"

    input:
    tuple val(sample), 
        val(ch_matrix), 
        val(job_id)
 
    output:
    tuple val(job_id), 
        val(sample), 
        path("afm.h5ad"), emit: afm
    
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
