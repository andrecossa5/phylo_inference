// INTERNAL_CONSISTENCY module
nextflow.enable.dsl = 2

process INTERNAL_CONSISTENCY {

    tag "${sample}_${filtering}_${boot_option}_${solver}_${metric}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering_key}/${solver}/${params.metric}/${boot_method}/", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering_key), 
        val(boot_input),
        val(boot_method),
        val(boot_replicates),
        val(solver),
        val(trees)
    
    output:
    tuple val(sample), 
        val(filtering_key), 
        val(boot_method),
        val(solver),
        path('trees.pickle'),
        path('extended_supports.csv'), emit: results

    script:
    """
    python ${baseDir}/bin/evaluate_tree/internal_evaluation.py \
    --sample_name ${sample} \
    --filtering ${filtering_key} \
    --solver ${solver} \
    --metric ${params.metric} \
    --boot_method ${boot_method} \
    --trees "${trees}"
    """

    stub:
    """
    touch trees.pickle
    touch extended_supports.csv
    """

}
