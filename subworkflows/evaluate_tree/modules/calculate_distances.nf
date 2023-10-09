// CALC_DISTANCES module
nextflow.enable.dsl = 2

// Module
process CALC_DISTANCES {

    tag "${sample}_${filtering}_${metric}"

    input: 
    tuple val(sample), 
        val(filtering), 
        val(metric),
        path(original_input)

    output:
    tuple val(sample), 
        val(filtering), 
        val(metric),
        path(original_input),
        path("dist.npz"), emit: dist
    
    script:
    """
    python ${baseDir}/bin/calculate_distances.py \
    -p ${original_input} \
    --metric ${metric} \
    --ncores ${task.cpus}
    """

    stub:
    """
    touch dist.npz
    """

}
