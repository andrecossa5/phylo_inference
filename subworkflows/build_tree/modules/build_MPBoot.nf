// BUILD_IQTREE module

nextflow.enable.dsl = 2

process BUILD_MPBOOT {

    tag "${sample}: ${filtering_key}, ${metric}, ${boot_method}, ${solver}, n=${boot_replicate}"

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
        path("rep_${boot_replicate}.newick"), emit: tree
    
    script:
    """
    mpboot -s ${input_folder}/sequences.fasta
    mv sequences.fasta.treefile rep_${boot_replicate}.newick
    """

    stub:
    """
    touch rep_${boot_replicate}.newick
    """

}