// BUILD_MPBOOT module

nextflow.enable.dsl = 2

process BUILD_MPBOOT {

    tag "${sample}: ${filtering_key}"

    input: 
    tuple val(sample), val(filtering_key), path(input_folder)

    output:
    tuple val(sample), val(filtering_key), path(input_folder), path("tree.newick"), emit: tree
    
    script:
    """
    mpboot -s ${input_folder}/sequences.fasta
    mv sequences.fasta.treefile tree.newick
    """

    stub:
    """
    touch tree.newick
    """

}
