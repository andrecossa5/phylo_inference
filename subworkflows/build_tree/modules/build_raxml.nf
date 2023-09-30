// BUILD_RAXML module

nextflow.enable.dsl = 2

process BUILD_RAXML {

    tag "${sample}_${filtering}_${boot_option}_${boot_replicate}_${solver}_${metric}"

    input:
    tuple val(sample), 
        val(filtering), 
        path(input_folder), 
        val(boot_replicate), 
        val(boot_option),
        path(bootstrapped_input),
        val(solver),
        val(metric)
    
    output:
    tuple val(sample), 
        val(filtering), 
        val(boot_replicate), 
        val(boot_option),
        val(solver),
        val(metric),
        path("rep${boot_replicate}.newick"), emit: tree
    
    script:
    """
    echo ${solver} > rep${boot_replicate}.txt
    """

    stub:
    """
    touch rep${boot_replicate}.txt
    """

}