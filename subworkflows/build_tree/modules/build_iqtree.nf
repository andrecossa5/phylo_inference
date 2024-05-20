// BUILD_IQTREE module

nextflow.enable.dsl = 2

process BUILD_IQTREE {

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
    iqtree -s ${input_folder}/sequences.fasta 
    python ${baseDir}/bin/build_tree/iqtree_to_newick.py ${boot_replicate} 
    """

    stub:
    """
    touch rep${boot_replicate}.newick
    """

}