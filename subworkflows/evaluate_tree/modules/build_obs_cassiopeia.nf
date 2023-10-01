// BUILD_OBSERVED_CASSIOPEIA module

nextflow.enable.dsl = 2


process BUILD_OBSERVED_CASSIOPEIA {

    tag "${sample}_${filtering}_${solver}_${metric}"

    input:
    tuple val(sample), 
        val(filtering), 
        path(input_folder), 
        val(solver),
        val(metric)
    
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
    --metric ${metric} \
    --solver ${solver} \
    --name tree \
    --ncores ${task.cpus}
    """

    stub:
    """
    touch tree.newick
    """

}
