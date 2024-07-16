// IQTREE module

nextflow.enable.dsl = 2

process IQTREE {

    tag "${sample}: ${filtering_key}"

    input: 
    tuple val(sample), val(filtering_key), path(input_folder)

    output:
    tuple val(sample), val(filtering_key), path(input_folder), path("final_tree.newick"), emit: tree
    
    script:
    """
    iqtree -s ${input_folder}/sequences.fasta 
    python ${baseDir}/bin/iqtree/iqtree_to_newick.py
    """

    stub:
    """
    touch final_tree.newick
    """

}