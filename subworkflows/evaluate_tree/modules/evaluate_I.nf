// EVALUATE_I module
nextflow.enable.dsl = 2

process EVALUATE_I {

    tag "${sample}_${filtering}_${boot_option}_${solver}_${metric}"

    // Publish
    publishDir "${params.outdir}/${sample}/${filtering}/${solver}/${metric}/${boot_option}/", mode: 'copy'

    input:
    tuple val(sample), 
        val(filtering), 
        val(solver),
        val(metric),
        val(boot_option),
        val(boot_trees),
        path(obs_tree)
    
    output:
    tuple val(sample), 
        val(filtering), 
        val(solver),
        val(metric),
        val(boot_option),
        path('trees.pickle'),
        path('extended_supports.csv'),
        path('results.txt'), emit: results
    
    script:
    """
    python ${baseDir}/bin/evaluation_I.py \
        --sample_name ${sample} \
        --filtering ${filtering} \
        --solver ${solver} \
        --metric ${metric} \
        --boot_method ${boot_option} \
        --obs_tree ${obs_tree} \
        --boot_trees "${boot_trees}" \
        --ncores ${task.cpus} 
    """

    stub:
    """
    touch results.txt
    touch trees.pickle
    touch extended_supports.csv
    """

}
