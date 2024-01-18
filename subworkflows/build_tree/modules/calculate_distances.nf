// CALC_DISTANCES module
nextflow.enable.dsl = 2

// Module
process CALC_DISTANCES {

    tag "${sample}_${filtering}_${metric}_${boot_option}_${boot_replicate}"

    input: 
    tuple val(sample), 
        val(filtering), 
        val(boot_option),
        val(boot_replicate), 
        val(metric),
        path(bootstrapped_input)

    output:
    tuple val(sample), 
        val(filtering), 
        val(boot_option),
        val(boot_replicate), 
        val(metric),
        path("dist.npz"), emit: dist
    
    script:
    """
    python ${baseDir}/bin/build_tree/calculate_distances.py \
    -p ${bootstrapped_input} \
    --metric ${metric} \
    --ncores ${task.cpus}
    """

    stub:
    """
    touch dist.npz
    """

}
