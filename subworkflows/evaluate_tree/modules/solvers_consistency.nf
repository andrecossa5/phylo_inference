//SOLVERS_CONSISTENCY module

nextflow.enable.dsl = 2

process SOLVERS_CONSISTENCY {

    tag "${sample}_${filtering}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering_key}/", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering_key), 
        val(boot_input),
        val(boot_method),
        val(boot_rep),
        val(solver),
        val(obs_tree)
    
    output: 
    tuple val(sample), 
        val(filtering_key), 
        path('results.csv'), emit: results
    
    script:
    """
    python ${baseDir}/bin/evaluate_tree/solvers_evaluation.py \
    --sample_name ${sample} \
    --filtering_key ${filtering_key} \
    --solver "${solver}" \
    --metric ${params.metric} \
    --obs_tree "${obs_tree}"
    """

    stub:
    """
    touch results.csv
    """

}

