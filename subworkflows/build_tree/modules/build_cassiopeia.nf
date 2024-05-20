// BUILD_CASSIOPEIA module

nextflow.enable.dsl = 2


process BUILD_CASSIOPEIA {

    tag "${sample}: ${filtering_key}, ${metric}, ${boot_method}, ${solver}, n=${boot_replicate}"

    input: 
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        val(boot_method),
        val(boot_replicate),
        val(solver),
        path(dist)

    output:
    tuple val(sample), 
        val(filtering_key), 
        path(input_folder),
        val(boot_method),
        val(boot_replicate),
        val(solver),
        path(dist),
        path("rep${boot_replicate}.newick"), emit: tree
    
    script:
    """
    python ${baseDir}/bin/build_tree/build_cassiopeia.py \
    -p ${input_folder} \
    -d ${dist} \
    --solver ${solver} \
    --name rep${boot_replicate} \
    --ncores ${task.cpus} 
    --path_filtering ${params.path_filtering} \
    --filtering_key ${filtering_key}
    """

    stub:
    """
    touch rep${boot_replicate}.newick
    """

}
