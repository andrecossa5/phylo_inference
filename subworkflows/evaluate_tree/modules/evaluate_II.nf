// EVALUATE_II module
nextflow.enable.dsl = 2

process EVALUATE_II {

    tag "${sample}_${filtering}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering}/", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering), 
        val(solver),
        val(metric),
        val(obs_tree)
    
    output:
    tuple val(sample), 
        val(filtering), 
        path('results.csv'), emit: results
    
    script:
    """
    python ${baseDir}/bin/evaluate_tree/evaluation_II.py \
    --sample_name ${sample} \
    --filtering ${filtering} \
    --solver "${solver}" \
    --metric "${metric}" \
    --obs_tree "${obs_tree}"
    """

    stub:
    """
    touch results.csv
    """

}

