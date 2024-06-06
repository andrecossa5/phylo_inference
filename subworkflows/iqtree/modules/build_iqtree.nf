// BUILD_IQTREE module

nextflow.enable.dsl = 2

process BUILD_IQTREE {

    tag "${sample}: ${filtering_key}"

    input: 
    tuple val(sample), val(filtering_key), path(input_folder)

    output:
    tuple val(sample), val(filtering_key), path(input_folder), path("tree.newick"), emit: tree
    
    script:
    """
    iqtree -s ${input_folder}/sequences.fasta 
    python ${baseDir}/bin/iqtree/iqtree_to_newick.py
    """

    stub:
    """
    touch tree.newick
    """

}