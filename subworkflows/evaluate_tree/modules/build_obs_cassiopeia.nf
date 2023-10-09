// BUILD_OBSERVED_CASSIOPEIA module

nextflow.enable.dsl = 2


process BUILD_OBSERVED_CASSIOPEIA {

    tag "${sample}_${filtering}_${solver}_${metric}"

    input:
    tuple val(sample), 
        val(filtering), 
        val(solver),
        val(metric),
        path(input_folder),
        path(dist)
    
    output:
    tuple val(sample), 
        val(filtering), 
        val(solver),
        val(metric),
        path("tree.newick"), emit: tree
    
    script:
    """
    python ${baseDir}/bin/build_cassiopeia.py \
    -p ${input_folder} \
    -d ${dist} \
    --solver ${solver} \
    --name tree \
    --ncores ${task.cpus} 
    """

    stub:
    """
    touch tree.newick
    """

}
