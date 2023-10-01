// BUILD_OBSERVED_SCITE module

nextflow.enable.dsl = 2


process BUILD_OBSERVED_SCITE {

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
    touch tree.newick
    """

    stub:
    """
    touch tree.newick
    """

}