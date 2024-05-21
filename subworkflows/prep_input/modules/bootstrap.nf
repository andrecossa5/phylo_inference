// BOOTSTRAP module

nextflow.enable.dsl = 2


process BOOTSTRAP {

    tag "${sample}: ${filtering_key}, ${boot_method}, rep ${boot_replicate}"

    input:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        val(boot_method),
        val(boot_replicate)
    
    output:
    tuple val(sample), 
        val(filtering_key),
        path(bootstrapped_input),
        val(boot_method),
        val(boot_replicate), emit: bootstrapped_input
    
    script:
    """
    python ${baseDir}/bin/prep_input/bootstrap.py \
    -p ${input_folder} \
    --method ${boot_method} \
    --feature_resampling_perc ${params.feature_resampling_perc} \
    --path_filtering ${params.path_filtering} \
    --filtering_key ${filtering_key} \
    --boot_replicate ${boot_replicate} 
    """

    stub:
    """
    mkdir bootstrapped_input
    """

}
