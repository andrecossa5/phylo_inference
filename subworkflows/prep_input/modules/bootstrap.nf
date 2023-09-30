// BOOTSTRAP module

nextflow.enable.dsl = 2


process BOOTSTRAP {

    tag "${sample}_${filtering}_${boot_option}_${boot_replicate}"

    input:
    tuple val(sample), val(filtering), path(input_folder), val(boot_replicate), val(boot_option)
    
    output:
    tuple val(sample), 
        val(filtering), 
        path(input_folder), 
        val(boot_replicate), 
        val(boot_option),
        path(bootstrapped_input), emit: bootstrapped_input
    
    script:
    """
    python ${baseDir}/bin/bootstrap.py \
        -p ${input_folder} \
        --method ${boot_option} 
    """

    stub:
    """
    mkdir bootstrapped_input
    """

}
