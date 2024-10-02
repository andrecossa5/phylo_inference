// IQTREE module

nextflow.enable.dsl = 2

process IQTREE {

    tag "${sample}: ${job_id}"

    input: 
    tuple val(job_id), val(sample), path(input_folder)

    output:
    tuple val(job_id), val(sample), path(input_folder), path("final_tree.newick"), emit: tree
    
    script:
    """
    iqtree -s ${input_folder}/genotypes.fa
    python ${baseDir}/build_tree/iqtree_to_newick.py
    """

    stub:
    """
    touch final_tree.newick
    """

}