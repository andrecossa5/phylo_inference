// MPBOOT module

nextflow.enable.dsl = 2

process MPBOOT {

    tag "${sample}: ${filtering_key}"

    input: 
    tuple val(sample), val(filtering_key), path(input_folder)

    output:
    tuple val(sample), val(filtering_key), path(input_folder), path("final_tree.newick"), emit: tree
    
    script:
    """
    mpboot -s ${input_folder}/sequences.fasta
    mv sequences.fasta.treefile final_tree.newick
    """

    stub:
    """
    touch final_tree.newick
    """

}
