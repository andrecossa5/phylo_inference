// BUILD_CASSIOPEIA module

nextflow.enable.dsl = 2


process BUILD_CASSIOPEIA {

    tag "${sample}_${filtering}_${boot_option}_${boot_replicate}_${solver}_${metric}"

    input: 
    tuple val(sample), 
        val(filtering), 
        val(boot_option),
        val(boot_replicate), 
        val(metric),
        val(solver),
        path(bootstrap_input),
        path(dist)

    output:
    tuple val(sample), 
        val(filtering), 
        val(boot_option),
        val(boot_replicate), 
        val(solver),
        val(metric),
        path("rep${boot_replicate}.newick"), emit: tree
    
    script:
    """
    python ${baseDir}/bin/build_tree/build_cassiopeia.py \
    -p ${bootstrap_input} \
    -d ${dist} \
    --solver ${solver} \
    --name rep${boot_replicate} \
    --ncores ${task.cpus} 
    """

    stub:
    """
    touch rep${boot_replicate}.newick
    """

}
