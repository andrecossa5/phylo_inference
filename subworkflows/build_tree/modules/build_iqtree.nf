// BUILD_IQTREE module

nextflow.enable.dsl = 2

process BUILD_IQTREE {

    tag "${sample}: ${filtering_key}, ${metric}, ${boot_method}, ${solver}, rep=${boot_replicate}"

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
    iqtree -s ${input_folder}/sequences.fasta 
    python ${baseDir}/bin/build_tree/iqtree_to_newick.py ${boot_replicate} 
    """

    stub:
    """
    touch rep_${boot_replicate}.newick
    """

}