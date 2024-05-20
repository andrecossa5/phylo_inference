// CALC_DISTANCES module
nextflow.enable.dsl = 2

// Module
process CALC_DISTANCES {

    tag "${sample}: ${filtering_key}, ${metric}, ${boot_method}, n=${boot_replicate}"

    input:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        val(boot_method),
        val(boot_replicate),
        val(solver)

    output:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        val(boot_method),
        val(boot_replicate),
        val(solver),
        path("dist.npz"), emit: dist
    
    script:
    """
    python ${baseDir}/bin/build_tree/calculate_distances.py \
    -p ${bootstrapped_input} \
    --metric ${metric} \
    --path_filtering ${params.path_filtering} \
    --filtering_key ${filtering_key}
    --ncores ${task.cpus}
    """

    stub:
    """
    touch dist.npz
    """

}
