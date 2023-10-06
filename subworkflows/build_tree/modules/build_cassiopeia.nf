// BUILD_CASSIOPEIA module

nextflow.enable.dsl = 2


process BUILD_CASSIOPEIA {

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
    python ${baseDir}/bin/build_cassiopeia.py \
    -p ${bootstrapped_input} \
    --metric ${metric} \
    --solver ${solver} \
    --name rep${boot_replicate} \
    --ncores ${task.cpus} 
    """

    stub:
    """
    touch rep${boot_replicate}.newick
    """

}
